
import os
import unittest
import numpy as np

from Lung_Segmentation import lung_segment


class LungSegmentationTest(unittest.TestCase):

    def setUp(self):
        test_data_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'resources', 'lung_image_132.dat')
        img_data = np.fromfile(test_data_file, 'int16')
        self._img_data = np.reshape(img_data, (307, 307))

    def test_lung_segment(self):
        seg, watershed, mark_int, mark_ext, mark_watershed, wat_left, wat_right = lung_segment(self._img_data)
        self.assertEqual((307, 307), mark_int.shape)
        self.assertFalse(mark_int[32, 48])
        self.assertEqual((307, 307), mark_ext.shape)
        self.assertEqual((307, 307), watershed.shape)


if __name__ == '__main__':
    unittest.main()
