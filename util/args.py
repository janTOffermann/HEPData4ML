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
        if isinstance(values, list):
            # Multiple arguments: -p 550 560 570
            try:
                result = [float(v) for v in values]
            except ValueError as e:
                parser.error("Invalid float value in {}: {}".format(option_string,e))
        else:
            # Single argument, might be comma/space separated: -p "550,560,570"
            if ',' in values:
                float_strings = [s.strip() for s in values.split(',')]
            else:
                float_strings = values.split()

            try:
                result = [float(s) for s in float_strings]
            except ValueError as e:
                parser.error("Invalid float value in {}: {}".format(option_string,e))

        setattr(namespace, self.dest, result)
