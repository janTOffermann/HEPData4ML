import argparse as ap

def none_or_str(value): # see https://stackoverflow.com/a/48295546
    if value == 'None':
        return None
    return value

def parse_mc_steps(value):
    """Parse MC steps from a list of strings, or a space- or comma-separated string"""
    ALLOWED_STEPS = ['generation', 'simulation', 'reconstruction']

    if isinstance(value, list):
        steps = value  # Already a list from nargs
    else:
        # Handle both comma and space separation
        if ',' in value:
            steps = [s.strip() for s in value.split(',')]
        else:
            steps = value.split()

    # Validate steps
    invalid = [step for step in steps if step not in ALLOWED_STEPS]
    if invalid:
        raise ap.ArgumentTypeError("Invalid steps: {}".format(invalid))
    return steps

def parse_float_list(value):
    """Parse a list of floats directly, or from a space- or comma-separated string"""

    if isinstance(value, list):
        floats = value  # Already a list from nargs
    else:
        # Handle both comma and space separation
        if ',' in value:
            float_strings = [s.strip() for s in value.split(',')]
        else:
            float_strings = value.split()

    # Convert to floats and validate
    try:
        floats = [float(s) for s in float_strings]
    except ValueError as e:
        raise ap.ArgumentTypeError("Invalid float value: {}".format(e))
    return floats