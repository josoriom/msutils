import unittest

from tests import support

MASS = 90.05550
FROM_RT = 0.0
TO_RT = 5.0
PPM_TOLERANCE = 20.0
MZ_TOLERANCE = 0.005

FILE_TOTAL = 34241
HEADER_BYTES = 1024

BYTES_TO_OPEN = FILE_TOTAL + HEADER_BYTES
REQUESTS_TO_OPEN = 2
BYTES_FOR_ONE_QUERY = 0
REQUESTS_FOR_ONE_QUERY = 0

BIG_FIXTURE = (
    support.REPO_ROOT / "core" / "tests" / "fixtures" / "check_features.ion"
)
BIG_REQUESTS_TO_OPEN = 2


@unittest.skipIf(support.library_problem, support.library_problem)
class PrefetchRanges(unittest.TestCase):
    def test_bytes_already_held_are_dropped_before_the_merge(self):
        source = support.recording_source(FILE_TOTAL)
        source.cache[(0, HEADER_BYTES)] = bytes(HEADER_BYTES)
        source.prefetch([(0, HEADER_BYTES), (2758, 64)])
        self.assertEqual(source.fetched, [(2758, 64)])

    def test_a_held_header_does_not_drag_in_the_gap_behind_it(self):
        source = support.recording_source(FILE_TOTAL)
        source.cache[(0, HEADER_BYTES)] = bytes(HEADER_BYTES)
        source.prefetch([(0, HEADER_BYTES), (2758, 64)])
        self.assertEqual(source.bytes_fetched, 64)

    def test_near_ranges_still_become_one_fetch(self):
        source = support.recording_source(FILE_TOTAL)
        source.prefetch([(100, 10), (140, 10)])
        self.assertEqual(source.fetched, [(100, 50)])

    def test_far_ranges_stay_apart(self):
        source = support.recording_source(FILE_TOTAL)
        source.prefetch([(100, 10), (300000, 10)])
        self.assertEqual(source.fetched, [(100, 10), (300000, 10)])

    def test_empty_ranges_are_ignored(self):
        source = support.recording_source(FILE_TOTAL)
        source.prefetch([(100, 0)])
        self.assertEqual(source.fetched, [])


@unittest.skipIf(support.library_problem, support.library_problem)
class SourceCache(unittest.TestCase):
    def test_a_missed_read_is_kept(self):
        source = support.recording_source(FILE_TOTAL)
        source.read(500, 20)
        source.read(500, 20)
        self.assertEqual(source.fetched, [(500, 20)])

    def test_a_read_inside_a_held_block_is_served_from_it(self):
        source = support.recording_source(FILE_TOTAL)
        source.read(500, 20)
        source.read(505, 5)
        self.assertEqual(source.fetched, [(500, 20)])


@unittest.skipIf(support.library_problem, support.library_problem)
class PlanSymbols(unittest.TestCase):
    def test_plan_open_is_bound(self):
        self.assertIn("plan_open", support._bridge._ABI._REQUIRED)
        self.assertTrue(hasattr(support.quantion._abi, "plan_open"))

    def test_plan_eic_is_bound(self):
        self.assertIn("plan_eic", support._bridge._ABI._REQUIRED)
        self.assertTrue(hasattr(support.quantion._abi, "plan_eic"))


@unittest.skipIf(support.library_problem, support.library_problem)
@unittest.skipIf(support.fixture_problem, support.fixture_problem)
class RemoteReading(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = support.RangeServer(support.ION_FIXTURE)

    @classmethod
    def tearDownClass(cls):
        cls.server.stop()

    def setUp(self):
        self.server.reset()

    def open_remote(self):
        remote = support.quantion.parse_ion_remote(self.server.url)
        self.addCleanup(remote.dispose)
        return remote

    def open_local(self):
        local = support.quantion.parse_ion_path(str(support.ION_FIXTURE))
        self.addCleanup(local.dispose)
        return local

    def eic_of(self, file):
        return support.quantion.calculate_eic(
            file, MASS, FROM_RT, TO_RT, PPM_TOLERANCE, MZ_TOLERANCE
        )

    def test_the_fixture_is_the_size_the_counts_were_measured_on(self):
        self.assertEqual(support.ION_FIXTURE.stat().st_size, FILE_TOTAL)

    def test_opening_fetches_the_planned_bytes_only(self):
        self.open_remote()
        self.assertEqual(self.server.bytes_sent, BYTES_TO_OPEN)

    def test_opening_takes_two_requests(self):
        self.open_remote()
        self.assertEqual(self.server.requests, REQUESTS_TO_OPEN)

    def test_opening_holds_the_whole_file(self):
        remote = self.open_remote()
        source = remote._remote_source
        self.assertIsNotNone(source.read_from_cache(0, FILE_TOTAL))

    def test_opening_reads_the_header_twice(self):
        remote = self.open_remote()
        source = remote._remote_source
        starts_at_zero = sorted(key for key in source.cache if key[0] == 0)
        self.assertEqual(starts_at_zero, [(0, HEADER_BYTES), (0, FILE_TOTAL)])

    def test_one_query_takes_no_request(self):
        remote = self.open_remote()
        self.server.reset()
        self.eic_of(remote)
        self.assertEqual(self.server.requests, REQUESTS_FOR_ONE_QUERY)

    def test_one_query_fetches_no_bytes(self):
        remote = self.open_remote()
        self.server.reset()
        self.eic_of(remote)
        self.assertEqual(self.server.bytes_sent, BYTES_FOR_ONE_QUERY)

    def test_a_repeated_query_asks_the_server_for_nothing(self):
        remote = self.open_remote()
        self.eic_of(remote)
        self.server.reset()
        self.eic_of(remote)
        self.assertEqual(self.server.requests, 0)

    def test_a_file_over_the_gap_is_opened_in_part(self):
        server = support.RangeServer(BIG_FIXTURE)
        self.addCleanup(server.stop)
        remote = support.quantion.parse_ion_remote(server.url)
        self.addCleanup(remote.dispose)
        self.assertEqual(server.requests, BIG_REQUESTS_TO_OPEN)
        self.assertLess(server.bytes_sent, BIG_FIXTURE.stat().st_size)

    def test_the_remote_answer_matches_the_local_one(self):
        remote_eic = self.eic_of(self.open_remote())
        local_eic = self.eic_of(self.open_local())
        self.assertEqual(remote_eic["x"].tolist(), local_eic["x"].tolist())
        self.assertEqual(remote_eic["y"].tolist(), local_eic["y"].tolist())


@unittest.skipIf(support.library_problem, support.library_problem)
@unittest.skipIf(support.fixture_problem, support.fixture_problem)
class ParseIonDispatch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = support.RangeServer(support.ION_FIXTURE)

    @classmethod
    def tearDownClass(cls):
        cls.server.stop()

    def open_and_keep(self, source):
        file = support.quantion.parse_ion(source)
        self.addCleanup(file.dispose)
        return file

    def test_a_url_goes_to_the_remote_reader(self):
        file = self.open_and_keep(self.server.url)
        self.assertTrue(hasattr(file, "_remote_source"))

    def test_a_path_goes_to_the_disk_reader(self):
        file = self.open_and_keep(str(support.ION_FIXTURE))
        self.assertFalse(hasattr(file, "_remote_source"))

    def test_bytes_go_to_the_memory_reader(self):
        file = self.open_and_keep(support.ION_FIXTURE.read_bytes())
        self.assertFalse(hasattr(file, "_remote_source"))

    def test_parse_ion_raw_refuses_a_path(self):
        with self.assertRaises(TypeError):
            support.quantion.parse_ion_raw(str(support.ION_FIXTURE))

    def test_parse_ion_remote_refuses_a_path(self):
        with self.assertRaises(ValueError):
            support.quantion.parse_ion_remote(str(support.ION_FIXTURE))
