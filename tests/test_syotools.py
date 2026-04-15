import os
import sys
import glob

import pytest
import numpy as np

from syotools.utils.yaml_utils import read_yaml, write_yaml
from syotools.wrappers.common import compute_observation, check_relative_diff


test_list = ["snrs", "exptimes", "magnitudes", "seds"]

'''
LOAD IT
'''

test_setups = []
test_names = []
for item in test_list:
    temp_files = glob.glob(f"{os.path.dirname(__file__)}/baselines/{item}/*.yml")
    for temp_file in temp_files:
        single_input = {"filename": temp_file}
        test_names.append(os.path.split(temp_file)[-1].replace(".yml", ""))
        test_setups.append(single_input)

def formatted_print(expected, result):
    """
    Formatting for the single-test runner version
    """
    print(f"                 Old:                 New: ")
    for idx, item in enumerate(expected):
        print(idx)
        for value in item:
            print(f" {value:>6}: {item[value]:>16}           {result[idx][value]:>16}")

@pytest.mark.parametrize("inputs", test_setups, ids = test_names)
def test_files(inputs):
    testfile = read_yaml(inputs["filename"])

    # write mode is the manual mode where you're trying to set (just tests without expected values) or reset (all tests)
    write = False
    if ("reset" in inputs and inputs["reset"]) or ("set" in inputs and inputs["set"]):
        write = True

    try:
        actual = compute_observation(testfile["telescope"], instrument=testfile["instrument"], sed=testfile["sed"], 
                    magnitude=testfile["magnitude"], snr=testfile["snr"], exptime=testfile["exptime"], 
                    redshift=testfile["redshift"], extinction=testfile["extinction"], target=testfile["target"])
    except Exception as err:
        if write: # this will never be true when run through pytest
            testfile["xfail"] = str(err)
            print("x", end="")
            write_yaml(testfile, inputs["filename"])
            return 1
        else:
            if "xfail" in testfile: # if this failed when first generated, it's an xfail
                if testfile["xfail"] == str(err): # Does it match that error?
                    pytest.xfail(testfile["xfail"])
                else: # if it doesn't match the other error, report it fresh
                    raise err 
            else:
                raise err # if not, it's a genuine error to be analyzed.

    result = []
    if actual is not None:
        for band in actual:
            result.append({"mean": np.nanmean(band.value), "median": np.nanmedian(band.value), "std": np.nanstd(band.value), "len": len(band)})
    # Two cases to set values:
    # 1. set is in the input command AND it's true AND there's no "expected" value in the file
    # 2. reset is in the input command AND it's true
    if write == False:
        if "expected" in testfile:
            # print mode is if you're manually running this code just to get the results of a single test.
            if "print" in inputs:
                formatted_print(testfile["expected"], result)
            assert check_relative_diff(result, testfile["expected"], 0.0005) #1e-3)
        else:
            assert False, f"No comparison in file {inputs['filename']}."
    elif ("reset" in inputs and inputs["reset"]) or ("set" in inputs and inputs["set"] and "expected" not in testfile): # only set if it doesn't have a record already
            testfile["expected"] = result
            print(".", end="")
            write_yaml(testfile, inputs["filename"])

if __name__ == "__main__":
    if len(sys.argv) > 0:
        additional = sys.argv[1]
        if additional in ["set", "reset"]:
            for test in test_setups:
                test[sys.argv[1]] = True
                test_files(test)
        else:
            test_names = [additional]
            test_files({"filename": additional, "print": True})
