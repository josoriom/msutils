import unittest

from tests import support

TARGET_RT = 1.0
RT_WINDOW = 0.2


@unittest.skipIf(support.library_problem, support.library_problem)
@unittest.skipIf(support.fixture_problem, support.fixture_problem)
class ChromatogramIndex(unittest.TestCase):
    def setUp(self):
        self.file = support.quantion.parse_ion_path(str(support.ION_FIXTURE))
        self.addCleanup(self.file.dispose)

    def peaks_for(self, item):
        return support.quantion.get_peaks_from_chrom(self.file, [item])

    def test_a_missing_index_is_refused(self):
        with self.assertRaisesRegex(ValueError, "idx"):
            self.peaks_for({"rt": TARGET_RT, "window": RT_WINDOW})

    def test_an_empty_index_is_refused(self):
        with self.assertRaisesRegex(ValueError, "idx"):
            self.peaks_for({"idx": None, "rt": TARGET_RT, "window": RT_WINDOW})

    def test_a_negative_index_is_refused(self):
        with self.assertRaisesRegex(ValueError, "idx"):
            self.peaks_for({"idx": -1, "rt": TARGET_RT, "window": RT_WINDOW})

    def test_the_message_names_the_item(self):
        with self.assertRaisesRegex(ValueError, "item 1"):
            support.quantion.get_peaks_from_chrom(
                self.file,
                [
                    {"idx": 0, "rt": TARGET_RT, "window": RT_WINDOW},
                    {"idx": -3, "rt": TARGET_RT, "window": RT_WINDOW},
                ],
            )

    def test_a_good_index_is_accepted(self):
        rows = self.peaks_for({"idx": 0, "rt": TARGET_RT, "window": RT_WINDOW})
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["index"], 0)

    def test_the_older_index_key_still_works(self):
        rows = self.peaks_for({"index": 0, "rt": TARGET_RT, "window": RT_WINDOW})
        self.assertEqual(rows[0]["index"], 0)

    def test_a_good_row_is_not_all_zeros(self):
        rows = self.peaks_for({"idx": 0, "rt": TARGET_RT, "window": RT_WINDOW})
        self.assertNotEqual(rows[0]["total_area"], 0.0)
