import unittest

from tests import support

MASS = 90.05550
FROM_RT = 0.0
TO_RT = 5.0
PPM_TOLERANCE = 20.0
MZ_TOLERANCE = 0.005

FILE_TOTAL = 34353
HEADER_BYTES = 1024
GAP_AFTER_HEADER = 1734

BYTES_TO_OPEN = 25680
REQUESTS_TO_OPEN = 3
BYTES_FOR_ONE_QUERY = 1734
REQUESTS_FOR_ONE_QUERY = 1


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
        source.prefetch([(100, 10), (30000, 10)])
        self.assertEqual(source.fetched, [(100, 10), (30000, 10)])

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

    def test_opening_takes_three_requests(self):
        self.open_remote()
        self.assertEqual(self.server.requests, REQUESTS_TO_OPEN)

    def test_opening_leaves_the_gap_behind_the_header_unread(self):
        remote = self.open_remote()
        source = remote._remote_source
        self.assertIsNone(source.read_from_cache(HEADER_BYTES, GAP_AFTER_HEADER))

    def test_opening_reads_the_header_once(self):
        remote = self.open_remote()
        source = remote._remote_source
        starts_at_zero = [key for key in source.cache if key[0] == 0]
        self.assertEqual(starts_at_zero, [(0, HEADER_BYTES)])

    def test_one_query_takes_one_request(self):
        remote = self.open_remote()
        self.server.reset()
        self.eic_of(remote)
        self.assertEqual(self.server.requests, REQUESTS_FOR_ONE_QUERY)

    def test_one_query_fetches_only_what_it_needs(self):
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

    def test_the_whole_file_is_never_fetched(self):
        remote = self.open_remote()
        self.eic_of(remote)
        self.assertLess(self.server.bytes_sent, FILE_TOTAL)

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
