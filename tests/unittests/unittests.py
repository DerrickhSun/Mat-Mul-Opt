from utils import *
from unittest import TestCase

"""
- For each operation, you should write tests to test  on matrices of different sizes.
- Keep in mind that the tests provided in the starter code are NOT comprehensive. That is, we strongly
advise you to modify them and add new tests.
- Hint: use dp_mc_matrix to generate dumbpy and numc matrices with the same data and use
      cmp_dp_nc_matrix to compare the results
"""
class TestAdd(TestCase):
    def test_small_add(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(100):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(i%9 + 1, i%9 + 1, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix(i%9 + 1, i%9 + 1, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
            self.assertTrue(is_correct)
            totaltime += speed_up
        for i in range(100):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(i%9 + 1, i%10 + 1, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix(i%9 + 1, i%10 + 1, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/200
        print_speedup(totaltime)
        try:
            nc.Matrix(3, 3) + nc.Matrix(2, 2)
            self.assertTrue(False)
        except ValueError as e:
            print(e)
        try:
            3 + nc.Matrix(2, 2)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            nc.Matrix(2, 2) + 3
            self.assertTrue(False)
        except TypeError as e:
            print(e)

    def test_medium_add(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(50):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%12 + 1) * 100, (i%12 + 1) * 100, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix((i%12 + 1) * 100, (i%12 + 1) * 100, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
            self.assertTrue(is_correct)
            totaltime += speed_up
        for i in range(50):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%11 + 1) * 100, (i%12 + 1) * 100, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix((i%11 + 1) * 100, (i%12 + 1) * 100, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/100
        print_speedup(totaltime)
        try:
            nc.Matrix(300, 300) + nc.Matrix(200, 200)
            self.assertTrue(False)
        except ValueError as e:
            print(e)
        try:
            3 + nc.Matrix(200, 200)
            self.assertTrue(False)
        except TypeError as e:
            print(e)

    def test_large_add(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(0):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%12 + 2) * 2000, (i%12 + 1) * 2000, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix((i%12 + 2) * 2000, (i%12 + 1) * 2000, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
            self.assertTrue(is_correct)
            totaltime += speed_up
        for i in range(0):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%11 + 2) * 2000, (i%12 + 2) * 2000, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix((i%11 + 2) * 2000, (i%12 + 2) * 2000, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/10
        print_speedup(totaltime)
        try:
            nc.Matrix(3000, 3000) + nc.Matrix(2000, 2000)
            self.assertTrue(False)
        except ValueError as e:
            print(e)
        try:
            3 + nc.Matrix(2000, 2000)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
            pass

class TestSub(TestCase):
    def test_small_sub(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(50):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(i%9 + 1, i%9 + 1, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix(i%9 + 1, i%9 + 1, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "sub")
            self.assertTrue(is_correct)
            totaltime += speed_up
        for i in range(50):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(i%9 + 1, i%10 + 1, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix(i%9 + 1, i%10 + 1, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "sub")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/100
        print_speedup(totaltime)
        try:
            nc.Matrix(3, 3) - nc.Matrix(2, 2)
            self.assertTrue(False)
        except ValueError as e:
            print(e)
        try:
            3 - nc.Matrix(2, 2)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            nc.Matrix(2, 2) - 2
            self.assertTrue(False)
        except TypeError as e:
            print(e)
            pass

    def test_medium_sub(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(50):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%12 + 1) * 100, (i%12 + 1) * 100, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix((i%12 + 1) * 100, (i%12 + 1) * 100, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "sub")
            self.assertTrue(is_correct)
            totaltime += speed_up
        for i in range(50):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%11 + 1) * 100, (i%12 + 1) * 100, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix((i%11 + 1) * 100, (i%12 + 1) * 100, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "sub")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/100
        print_speedup(totaltime)
        try:
            nc.Matrix(300, 300) - nc.Matrix(200, 200)
            self.assertTrue(False)
        except ValueError as e:
            print(e)
        try:
            30 - nc.Matrix(200, 200)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            nc.Matrix(200, 200)-300
            self.assertTrue(False)
        except TypeError as e:
            print(e)
            pass

    def test_large_sub(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(0):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%12 + 1) * 1000, (i%12 + 1) * 1000, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix((i%12 + 1) * 1000, (i%12 + 1) * 1000, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
            self.assertTrue(is_correct)
            totaltime += speed_up
        for i in range(0):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%4 + 2) * 2000, (i%12 + 1) * 2000, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix((i%4 + 2) * 2000, (i%12 + 1) * 2000, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "add")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/10
        print_speedup(totaltime)
        try:
            nc.Matrix(3000, 3000) - nc.Matrix(2000, 2000)
            self.assertTrue(False)
        except ValueError as e:
            print(e)
        try:
            3 - nc.Matrix(2000, 2000)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
            pass

class TestAbs(TestCase):
    def test_small_abs(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(50):
            dp_mat, nc_mat = rand_dp_nc_matrix(i%9 + 1, i%9 + 1, seed=i)
            is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
            self.assertTrue(is_correct)
            totaltime += speed_up
        for i in range(50):
            dp_mat, nc_mat = rand_dp_nc_matrix(i%9 + 1, i%10 + 1, seed=i)
            is_correct, speed_up = compute([dp_mat], [nc_mat], "abs")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/100
        print_speedup(totaltime)
        pass

    def test_medium_abs(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(20):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%12 + 1) * 100, (i%12 + 1) * 100, seed=i)
            is_correct, speed_up = compute([dp_mat1], [nc_mat1], "abs")
            self.assertTrue(is_correct)
            totaltime += speed_up
        for i in range(20):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%11 + 1) * 100, (i%12 + 1) * 100, seed=i)
            is_correct, speed_up = compute([dp_mat1], [nc_mat1], "abs")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/40
        print_speedup(totaltime)
        pass

    def test_large_abs(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(0):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%12 + 2) * 2000, (i%12 + 2) * 2000, seed=i)
            is_correct, speed_up = compute([dp_mat1], [nc_mat1], "abs")
            self.assertTrue(is_correct)
            totaltime += speed_up
        for i in range(0):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%4 + 2) * 2000, (i%12 + 2) * 2000, seed=i)
            is_correct, speed_up = compute([dp_mat1], [nc_mat1], "abs")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/10
        print_speedup(totaltime)
        pass

class TestNeg(TestCase):
    def test_small_neg(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(50):
            dp_mat, nc_mat = rand_dp_nc_matrix(i%9 + 1, i%9 + 1, seed=i)
            is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
            self.assertTrue(is_correct)
            totaltime += speed_up
        for i in range(50):
            dp_mat, nc_mat = rand_dp_nc_matrix(i%9 + 1, i%10 + 1, seed=i)
            is_correct, speed_up = compute([dp_mat], [nc_mat], "neg")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/100
        print_speedup(totaltime)
        pass
    def test_medium_neg(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(20):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%12 + 1) * 100, (i%12 + 1) * 100, seed=i)
            is_correct, speed_up = compute([dp_mat1], [nc_mat1], "neg")
            self.assertTrue(is_correct)
            totaltime += speed_up
        for i in range(20):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%11 + 1) * 100, (i%12 + 1) * 100, seed=i)
            is_correct, speed_up = compute([dp_mat1], [nc_mat1], "neg")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/40
        print_speedup(totaltime)
        pass

    def test_large_neg(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(0):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%11 + 2) * 2000, (i%11 + 2) * 2000, seed=i)
            is_correct, speed_up = compute([dp_mat1], [nc_mat1], "neg")
            self.assertTrue(is_correct)
            totaltime += speed_up
        for i in range(0):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%11 + 1) * 2000, (i%12 + 1) * 2000, seed=i)
            is_correct, speed_up = compute([dp_mat1], [nc_mat1], "neg")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/10
        print_speedup(totaltime)
        pass

class TestMul(TestCase):
    def test_small_mul(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(100):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(i%9 + 1, i%9 + 1, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix(i%9 + 1, i%9 + 1, seed=i)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
            if not is_correct:
                print(i)
            self.assertTrue(is_correct)
            totaltime += speed_up
        for i in range(100):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(i%9 + 1, i%10 + 1, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix(i%10 + 1, i%11 + 1, seed=i)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
            if not is_correct:
                print(i)
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/200
        print_speedup(totaltime)
        try:
            nc.Matrix(3, 3) * nc.Matrix(2, 2)
            self.assertTrue(False)
        except ValueError as e:
            print(e)
        try:
            3 * nc.Matrix(2, 2)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            nc.Matrix(2, 2) * 3
            self.assertTrue(False)
        except TypeError as e:
            print(e)
            pass

    def test_medium_mul(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(10):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%12 + 1) * 100, (i%12 + 1) * 100, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix((i%12 + 1) * 100, (i%12 + 1) * 100, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
            self.assertTrue(is_correct)
            totaltime += speed_up
        for i in range(0):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%11 + 1) * 100, (i%12 + 1) * 100, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix((i%11 + 1) * 100, (i%12 + 1) * 100, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/20
        print_speedup(totaltime)
        try:
            nc.Matrix(30, 30) * nc.Matrix(20, 20)
            self.assertTrue(False)
        except ValueError as e:
            print(e)
        try:
            30 * nc.Matrix(20, 20)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            nc.Matrix(20, 20) * 30
            self.assertTrue(False)
        except TypeError as e:
            print(e)
            pass

    def test_large_mul(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(0):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%12 + 1) * 2000, (i%12 + 2) * 2000, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix((i%12 + 2) * 2000, (i%12 + 1) * 2000, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
            self.assertTrue(is_correct)
            totaltime += speed_up
        for i in range(0):
            dp_mat1, nc_mat1 = rand_dp_nc_matrix((i%11 + 2) * 2000, (i%12 + 2) * 2000, seed=i)
            dp_mat2, nc_mat2 = rand_dp_nc_matrix((i%11 + 2) * 2000, (i%12 + 2) * 2000, seed=i+1)
            is_correct, speed_up = compute([dp_mat1, dp_mat2], [nc_mat1, nc_mat2], "mul")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/10
        print_speedup(totaltime)
        pass

class TestPow(TestCase):
    def test_small_pow(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(50):
            dp_mat, nc_mat = rand_dp_nc_matrix(i%9 + 1, i%9 + 1, seed=i)
            is_correct, speed_up = compute([dp_mat, i%20], [nc_mat, i%20], "pow")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/500
        print_speedup(totaltime)
        try:
            nc.Matrix(3, 3) ** nc.Matrix(2, 2)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            3 ** nc.Matrix(2, 2)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            nc.Matrix(2, 2) ** -3
            self.assertTrue(False)
        except ValueError as e:
            print(e)
        try:
            nc.Matrix(2, 1) ** 3
            self.assertTrue(False)
        except ValueError as e:
            print(e)
            pass

    def test_medium_pow(self):
        # TODO: YOUR CODE HERE
        totaltime = 0
        for i in range(10):
            dp_mat, nc_mat = rand_dp_nc_matrix((i%12 + 1) * 100, (i%12 + 1) * 100, seed=i)
            is_correct, speed_up = compute([dp_mat, i%50], [nc_mat, i%50], "pow")
            self.assertTrue(is_correct)
            totaltime += speed_up
        
        totaltime = totaltime/100
        print_speedup(totaltime)
        try:
            nc.Matrix(30, 30) ** nc.Matrix(2, 2)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            30 ** nc.Matrix(20, 20)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            nc.Matrix(20, 20) ** -30
            self.assertTrue(False)
        except ValueError as e:
            print(e)
        try:
            nc.Matrix(20, 10) ** 30
            self.assertTrue(False)
        except ValueError as e:
            print(e)
            pass

    def test_large_pow(self):
        # TODO: YOUR CODE HERE
        pass

class TestGet(TestCase):
    def test_get(self):
        # TODO: YOUR CODE HERE
        for i in range(20): 
            dp_mat, nc_mat = rand_dp_nc_matrix(i%17 + 1, i%16 + 1, seed=i)
            rand_row = np.random.randint(dp_mat.shape[0])
            rand_col = np.random.randint(dp_mat.shape[1])
            self.assertEqual(round(dp_mat.get(rand_row, rand_col), decimal_places),
                round(nc_mat.get(rand_row, rand_col), decimal_places))
        try:
            nc_mat.get(1, 2, 3)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            nc_mat.get(1.0, 0)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            nc_mat.get(9, "hi")
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            nc_mat.get(9)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            nc_mat.get(-1, 0)
            self.assertTrue(False)
        except IndexError as e:
            print(e)
        try:
            nc_mat.get(9999999, 999999)
            self.assertTrue(False)
        except IndexError as e:
            print(e)
            pass

class TestSet(TestCase):
    def test_set(self):
        # TODO: YOUR CODE HERE
        
        for i in range(50):
            dp_mat, nc_mat = rand_dp_nc_matrix(i%17 + 1, i%16 + 1, seed=i)
            rand_row = np.random.randint(dp_mat.shape[0])
            rand_col = np.random.randint(dp_mat.shape[1])
            dp_mat.set(rand_row, rand_col, 2)
            nc_mat.set(rand_row, rand_col, 2)
            self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
        try:
            nc_mat.get(0, 0, 3, 0)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            nc_mat.get(0.0, 0, 1)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            nc_mat.get(0, "hi", 2)
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            nc_mat.get(0, 0, "wut")
            self.assertTrue(False)
        except TypeError as e:
            print(e)
        try:
            nc_mat.get(-1, 0, 9.0)
            self.assertTrue(False)
        except IndexError as e:
            print(e)
        try:
            nc_mat.get(9999999, 999999, 1)
            self.assertTrue(False)
        except IndexError as e:
            print(e)
            pass

class TestShape(TestCase):
    def test_shape(self):
        # TODO: YOUR CODE HERE
        for i in range(20):
            dp_mat, nc_mat = rand_dp_nc_matrix(i+2, i + 1, seed=i)
            self.assertTrue(dp_mat.shape == nc_mat.shape)

class TestIndexGet(TestCase):
    def test_index_get(self):
        # TODO: YOUR CODE HERE
        for i in range(20):
            dp_mat, nc_mat = rand_dp_nc_matrix(i+2, i + 1, seed=i)
            rand_row = np.random.randint(dp_mat.shape[0])
            rand_col = np.random.randint(dp_mat.shape[1])
            self.assertEqual(round(dp_mat[rand_row][rand_col], decimal_places),
                round(nc_mat[rand_row][rand_col], decimal_places))

class TestIndexSet(TestCase):
    def test_set(self):
        # TODO: YOUR CODE HERE
        for i in range(20):
            dp_mat, nc_mat = rand_dp_nc_matrix(i+2, i+2, seed=0)
            rand_row = np.random.randint(dp_mat.shape[0])
            rand_col = np.random.randint(dp_mat.shape[1])
            dp_mat[rand_row][rand_col] = 2
            nc_mat[rand_row][rand_col] = 2
            self.assertTrue(cmp_dp_nc_matrix(dp_mat, nc_mat))
            self.assertEquals(nc_mat[rand_row][rand_col], 2)

class TestSlice(TestCase):
    def test_slice(self):
        # TODO: YOUR CODE HERE
        for i in range(20):
            dp_mat, nc_mat = rand_dp_nc_matrix(i + 2, i + 2, seed=i)
            dp_mat1, nc_mat1 = rand_dp_nc_matrix(i + 2, i + 2, seed=i)
            self.assertTrue(cmp_dp_nc_matrix(dp_mat[0], nc_mat[0]))
            self.assertTrue(cmp_dp_nc_matrix(dp_mat[1], nc_mat[1]))
