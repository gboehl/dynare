function failed_tests = my_assert(failed_tests, success, test_name)
if ~success
    failed_tests = cat(1, test_failed, test_name);
end