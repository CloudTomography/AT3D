import copy
import numpy as np
import pyshdom
from unittest import TestCase

class SortTest(TestCase):
    @classmethod
    def setUpClass(cls):
        np.random.seed(3)

        x = np.random.uniform(low=0.0,high=100.0, size=2000).astype(np.int)
        y = np.random.uniform(low=0.0,high=100.0, size=2000).astype(np.float32)
        yi = np.random.uniform(low=0.0, high=100.0, size=2000).astype(np.int)

        x_sorted,y_sorted, yi_sorted = pyshdom.core.quicksort_new(x=copy.deepcopy(x),y=copy.deepcopy(y),
                                                       yi=copy.deepcopy(yi),
                                                       n=x.size,kflag=2)
        assert np.all(x[np.argsort(x)]== x_sorted)
        tests = []
        for i in range(x.min(),x.max()+1):
            tests.append(np.allclose(np.sort(y[np.where(x==i)]),
                        np.sort(y_sorted[np.where(x_sorted==i)])))
        assert np.all(tests)
        tests2 = []
        for i in range(x.min(),x.max()+1):
            tests2.append(np.allclose(np.sort(yi[np.where(x==i)]),
                        np.sort(yi_sorted[np.where(x_sorted==i)])))
        assert np.all(tests2)
        cls.x = x
        cls.y = y
        cls.yi = yi
        cls.x_sorted = x_sorted
        cls.y_sorted = y_sorted
        cls.yi_sorted = yi_sorted
    def test_x(self):
        self.assertTrue(np.all(self.x[np.argsort(self.x)]== self.x_sorted))
    def test_y(self):
        tests = []
        for i in range(self.x.min(),self.x.max()+1):
            tests.append(np.allclose(np.sort(self.y[np.where(self.x==i)]),
                        np.sort(self.y_sorted[np.where(self.x_sorted==i)])))
        self.assertTrue(np.all(tests))
    def test_yi(self):
        tests = []
        for i in range(self.x.min(),self.x.max()+1):
            tests.append(np.allclose(np.sort(self.yi[np.where(self.x==i)]),
                        np.sort(self.yi_sorted[np.where(self.x_sorted==i)])))
        self.assertTrue(tests)
