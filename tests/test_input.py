from insilico import validate_input


def test_valid_input():
    """
    Test that the input is properly check.
    """
    file_path = "tests/tests_files/input_test.yml"
    inp = validate_input(file_path)

    assert isinstance(inp, dict)
