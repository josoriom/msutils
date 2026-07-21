import ctypes
import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock

from tests import support

OLDER_VERSION = "0.9.0"
NEWER_VERSION = "0.10.0"
FUTURE_VERSION = "9.9.9"


@unittest.skipIf(support.library_problem, support.library_problem)
class LibraryLookupOrder(unittest.TestCase):
    def setUp(self):
        folder = TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        chosen_root = mock.patch.dict(
            os.environ, {"QUANTION_ARTIFACTS_ROOT": str(self.root)}
        )
        chosen_root.start()
        self.addCleanup(chosen_root.stop)
        self.own_version = support._bridge._own_version()

    def found(self):
        return support._bridge._bundled_lib_path()

    def test_an_unversioned_library_wins(self):
        wanted = support.place_library(self.root)
        support.place_library(self.root, self.own_version)
        support.place_library(self.root, NEWER_VERSION)
        support.place_manifest(self.root, NEWER_VERSION)
        self.assertEqual(self.found(), wanted)

    def test_an_unversioned_library_needs_no_manifest(self):
        wanted = support.place_library(self.root)
        self.assertEqual(self.found(), wanted)

    def test_the_own_version_needs_no_manifest(self):
        wanted = support.place_library(self.root, self.own_version)
        self.assertEqual(self.found(), wanted)

    def test_the_own_version_beats_a_higher_version(self):
        wanted = support.place_library(self.root, self.own_version)
        support.place_library(self.root, FUTURE_VERSION)
        support.place_manifest(self.root, FUTURE_VERSION)
        self.assertEqual(self.found(), wanted)

    def test_the_scan_reads_versions_as_numbers(self):
        support.place_library(self.root, OLDER_VERSION)
        support.place_manifest(self.root, OLDER_VERSION)
        wanted = support.place_library(self.root, NEWER_VERSION)
        support.place_manifest(self.root, NEWER_VERSION)
        self.assertEqual(self.found(), wanted)

    def test_the_scan_skips_a_version_without_a_manifest(self):
        support.place_library(self.root, FUTURE_VERSION)
        wanted = support.place_library(self.root, NEWER_VERSION)
        support.place_manifest(self.root, NEWER_VERSION)
        self.assertEqual(self.found(), wanted)

    def test_the_scan_skips_a_version_without_this_platform(self):
        support.place_manifest(self.root, FUTURE_VERSION)
        wanted = support.place_library(self.root, NEWER_VERSION)
        support.place_manifest(self.root, NEWER_VERSION)
        self.assertEqual(self.found(), wanted)

    def test_an_empty_root_falls_through_to_the_packaged_library(self):
        packaged = self.root / "packaged" / support._bridge._default_lib_name()
        with mock.patch.object(
            support._bridge, "_packaged_lib_path", lambda: packaged
        ):
            self.assertEqual(self.found(), packaged)

    def test_a_missing_root_falls_through_to_the_packaged_library(self):
        missing = self.root / "gone"
        with mock.patch.dict(os.environ, {"QUANTION_ARTIFACTS_ROOT": str(missing)}):
            packaged = self.root / "packaged" / support._bridge._default_lib_name()
            with mock.patch.object(
                support._bridge, "_packaged_lib_path", lambda: packaged
            ):
                self.assertEqual(self.found(), packaged)

    def test_a_root_with_only_a_manifest_falls_through(self):
        support.place_manifest(self.root, NEWER_VERSION)
        packaged = self.root / "packaged" / support._bridge._default_lib_name()
        with mock.patch.object(
            support._bridge, "_packaged_lib_path", lambda: packaged
        ):
            self.assertEqual(self.found(), packaged)


@unittest.skipIf(support.library_problem, support.library_problem)
class LibraryFromTheEnvironment(unittest.TestCase):
    def missing_path(self) -> str:
        folder = TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        return str(Path(folder.name) / support._bridge._default_lib_name())

    def test_the_env_path_is_tried_before_the_artifacts(self):
        asked_for = self.missing_path()
        chosen = {
            "QUANTION_LIB": asked_for,
            "QUANTION_ARTIFACTS_ROOT": str(support.ARTIFACTS_ROOT),
        }
        with mock.patch.dict(os.environ, chosen):
            with self.assertRaises(OSError) as caught:
                support._bridge.load_library()
        self.assertIn(asked_for, str(caught.exception))

    def test_a_given_path_is_tried_before_the_environment(self):
        asked_for = self.missing_path()
        with mock.patch.dict(os.environ, {"QUANTION_LIB": "/never/used"}):
            with self.assertRaises(OSError) as caught:
                support._bridge.load_library(asked_for)
        self.assertIn(asked_for, str(caught.exception))

    def test_the_old_env_name_is_not_read(self):
        empty_root = TemporaryDirectory()
        self.addCleanup(empty_root.cleanup)
        chosen = {
            "MSUTILS_LIB": self.missing_path(),
            "QUANTION_ARTIFACTS_ROOT": empty_root.name,
        }
        with mock.patch.dict(os.environ, chosen):
            os.environ.pop("QUANTION_LIB", None)
            with self.assertRaises(FileNotFoundError) as caught:
                support._bridge.load_library()
        self.assertIn("QUANTION_LIB", str(caught.exception))

    def test_the_missing_library_message_names_every_escape_hatch(self):
        empty_root = TemporaryDirectory()
        self.addCleanup(empty_root.cleanup)
        with mock.patch.dict(os.environ, {"QUANTION_ARTIFACTS_ROOT": empty_root.name}):
            os.environ.pop("QUANTION_LIB", None)
            with mock.patch.object(support._bridge, "_packaged_lib_path", lambda: None):
                with self.assertRaises(FileNotFoundError) as caught:
                    support._bridge.load_library()
        message = str(caught.exception)
        self.assertIn("QUANTION_LIB", message)
        self.assertIn("QUANTION_ARTIFACTS_ROOT", message)


@unittest.skipIf(support.library_problem, support.library_problem)
class LibraryHandshake(unittest.TestCase):
    def test_the_abi_version_symbol_is_the_renamed_one(self):
        self.assertEqual(
            support.quantion._abi.quantion_abi_version(),
            support._bridge.QUANTION_ABI_VERSION,
        )

    def test_the_struct_size_symbol_is_the_renamed_one(self):
        self.assertEqual(
            support.quantion._abi.quantion_sizeof_peak_options(),
            ctypes.sizeof(support._bridge.PeakOptions),
        )

    def test_the_old_symbol_names_are_gone(self):
        self.assertNotIn("msutils_abi_version", support._bridge._ABI._REQUIRED)
        self.assertNotIn(
            "msutils_sizeof_peak_options", support._bridge._ABI._REQUIRED
        )
        self.assertFalse(hasattr(support._bridge, "MSUTILS_ABI_VERSION"))

    def test_the_library_name_carries_the_new_name(self):
        self.assertIn("quantion", support._bridge._default_lib_name())


@unittest.skipIf(support.library_problem, support.library_problem)
class PackageRootDetection(unittest.TestCase):
    def setUp(self):
        folder = TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)

    def write_manifest(self, text: str) -> None:
        (self.root / "pyproject.toml").write_text(text, encoding="utf-8")

    def test_a_quantion_manifest_marks_the_root(self):
        self.write_manifest('[project]\nname = "quantion"\n')
        self.assertTrue(support._bridge._is_package_root(self.root))

    def test_another_manifest_does_not_mark_the_root(self):
        self.write_manifest('[project]\nname = "something-else"\n')
        self.assertFalse(support._bridge._is_package_root(self.root))

    def test_no_manifest_does_not_mark_the_root(self):
        self.assertFalse(support._bridge._is_package_root(self.root))
