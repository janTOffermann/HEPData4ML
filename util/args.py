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

class FloatListAction(ap.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if not isinstance(values, list):
            # This shouldn't happen with nargs='*', but handle it just in case
            values = [values]

        result = []
        for value in values:
            if isinstance(value, (int, float)):
                # Already a number
                result.append(float(value))
            elif isinstance(value, str):
                # String that might contain delimiters
                if ',' in value:
                    # Comma-separated
                    float_strings = [s.strip() for s in value.split(',')]
                else:
                    # Space-separated
                    float_strings = value.split()

                try:
                    result.extend([float(s) for s in float_strings if s.strip()])
                except ValueError as e:
                    parser.error("Invalid float value in {}: {}".format(option_string, e))
            else:
                parser.error("Unexpected value type in {}: {}".format(option_string, type(value)))

        setattr(namespace, self.dest, result)

