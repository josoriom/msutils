import inspect
import unittest
import warnings

from tests import support


@unittest.skipIf(support.library_problem, support.library_problem)
class RenamedClass(unittest.TestCase):
    def test_the_new_name_is_a_class(self):
        self.assertTrue(inspect.isclass(support.quantion.Quantion))

    def test_the_new_name_is_advertised(self):
        self.assertIn("Quantion", support.quantion.__all__)

    def test_the_old_name_still_resolves(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            old_name = support.quantion.MsUtils
        self.assertIs(old_name, support.quantion.Quantion)

    def test_the_old_name_warns(self):
        with warnings.catch_warnings(record=True) as seen:
            warnings.simplefilter("always")
            support.quantion.MsUtils
        self.assertEqual(len(seen), 1)
        self.assertTrue(issubclass(seen[0].category, DeprecationWarning))

    def test_the_warning_names_both_names(self):
        with warnings.catch_warnings(record=True) as seen:
            warnings.simplefilter("always")
            support.quantion.MsUtils
        message = str(seen[0].message)
        self.assertIn("MsUtils", message)
        self.assertIn("Quantion", message)

    def test_the_new_name_does_not_warn(self):
        with warnings.catch_warnings(record=True) as seen:
            warnings.simplefilter("always")
            support.quantion.Quantion
        self.assertEqual(seen, [])

    def test_the_old_name_is_not_advertised(self):
        self.assertNotIn("MsUtils", support.quantion.__all__)

    def test_an_unknown_name_still_fails(self):
        with self.assertRaises(AttributeError):
            support.quantion.no_such_name


@unittest.skipIf(support.library_problem, support.library_problem)
class RemovedStderrCapture(unittest.TestCase):
    def test_the_bridge_has_no_capture_pair(self):
        self.assertFalse(hasattr(support._bridge, "start_stderr_capture"))
        self.assertFalse(hasattr(support._bridge, "stop_stderr_capture"))

    def test_the_package_has_no_capture_pair(self):
        self.assertFalse(hasattr(support.quantion, "start_stderr_capture"))
        self.assertFalse(hasattr(support.quantion, "stop_stderr_capture"))

    def test_the_capture_pair_is_not_advertised(self):
        self.assertNotIn("start_stderr_capture", support.quantion.__all__)
        self.assertNotIn("stop_stderr_capture", support.quantion.__all__)


@unittest.skipIf(support.library_problem, support.library_problem)
class PublicSurface(unittest.TestCase):
    def test_every_advertised_name_exists(self):
        for name in support.quantion.__all__:
            self.assertTrue(hasattr(support.quantion, name), name)

    def test_the_remote_reader_replaces_the_url_reader(self):
        self.assertTrue(hasattr(support.quantion, "parse_ion_remote"))
        self.assertFalse(hasattr(support.quantion, "parse_ion_url"))

    def test_the_class_offers_the_remote_reader(self):
        self.assertTrue(hasattr(support.quantion.Quantion, "parse_ion_remote"))
        self.assertTrue(hasattr(support.quantion.Quantion, "parse_ion_raw"))
