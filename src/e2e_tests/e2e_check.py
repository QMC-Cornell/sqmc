"""Extracts energy from the output file and check against a reference file."""
import sys
import re

DETERMINISTIC_TOLERANCE = 0.01
STOCHASTIC_TOLERANCE = 0.05
COLORS = {"OKGREEN": "\033[92m", "WARNING": "\033[93m", "ENDC": "\033[0m"}

CHECKS = []
CHECKS.append({
    "title": "Variational Energy",
    "pattern": "Variational energy.*=\s*([-+]?[0-9]*\.?[0-9]+)",
    "tolerance": DETERMINISTIC_TOLERANCE
})
CHECKS.append({
    "title": "PT Correction",
    "pattern": "Second-order PT energy lowering.*=\s*([-+]?[0-9]*\.?[0-9]+)",
    "tolerance": DETERMINISTIC_TOLERANCE
})
CHECKS.append({
    "title": "PT Correction Uncertainty",
    "pattern": "Second-order PT energy lowering.*=\s*[-+]?[0-9]*\.?[0-9]+\s\+\-\s([0-9]*\.?[0-9]+)",
    "tolerance": STOCHASTIC_TOLERANCE
})

MATCH_OUTPUT = COLORS['OKGREEN'] + 'Match.' + COLORS['ENDC']
UNMATCH_OUTPUT = COLORS['WARNING'] + 'Not match.' + COLORS['ENDC']

def main():
    # Obtain output and reference filenames.
    argv = sys.argv
    if len(argv) != 3:
        print "Usage: python e2e_check.py [output] [reference]"
        return
    output_filename = argv[1]
    reference_filename = argv[2]

    output_values = {}
    reference_values = {}

    for check in CHECKS:
        print "%s:" % check['title']
        pattern = re.compile(check['pattern'])

        output_value_found = False
        reference_value_found = False

        for i, line in enumerate(open(output_filename)):
            for match in re.finditer(pattern, line):
                print "Found in %s line %s: %s" % (output_filename, i + 1, match.group(0))
                output_value = float(match.group(1))
                output_value_found = True
            if output_value_found:
                break

        for i, line in enumerate(open(reference_filename)):
            for match in re.finditer(pattern, line):
                print "Found in %s line %s: %s" % (reference_filename, i + 1, match.group(0))
                reference_value = float(match.group(1))
                reference_value_found = True
            if reference_value_found:
                break

        if output_value_found and reference_value_found:
            diff = abs(output_value - reference_value)
            relative_diff = abs(diff * 1.0 / output_value)
            if relative_diff < check['tolerance']:
                status = MATCH_OUTPUT
            else:
                status = UNMATCH_OUTPUT
            print "Diff: %f (%f %%) %s" % (diff, relative_diff * 100, status)
        elif not(output_value_found) and not(reference_value_found):
            print "Both not found. %s" % (MATCH_OUTPUT,)
        else:
            print UNMATCH_OUTPUT

        if output_value_found:
            output_values[check['title']] = output_value
        if reference_value_found:
            reference_values[check['title']] = reference_value

        print "-" * 80

    if 'PT Correction' in output_values \
        and 'PT Correction' in reference_values \
        and 'PT Correction Uncertainty' in reference_values:
        diff = abs(output_values['PT Correction'] - reference_values['PT Correction'])
        uncertainty = reference_values['PT Correction Uncertainty']
        print 'PT Correction is %f sigma away from reference value.' % (diff / uncertainty)
        print "-" * 80

if __name__ == "__main__":
    main()



