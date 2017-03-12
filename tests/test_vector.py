import setup_miso
import pysplicing

def test_vector_median():
    v1 = tuple(range(0, 9))
    assert(pysplicing.i_test_vector_median(v1) == 4.0)
    v2 = tuple(range(0, 10))
    assert(pysplicing.i_test_vector_median(v2) == 4.5)
    v3 = (1,)
    assert(pysplicing.i_test_vector_median(v3) == 1.0)
    v3 = (1,2)
    assert(pysplicing.i_test_vector_median(v3) == 1.5)    
    v4 = tuple(range(8, -1, -1))
    assert(pysplicing.i_test_vector_median(v4) == 4.0)
