"""
warnings.py

All the warnings
"""

# FIXME: it would be better to have some structure to he warning codes

WARNINGS = {
    "W001":
        "Record name: {} is not very descriptive",
    "W002":
        "Record has no product type",
    "W003":
        '"{}" is not a valid product type. Options are: {}',
    "W004":
        "No api value provided",
    "W005":
        "API value: {api} seems unlikely",
    "W011":
        "{} date format: {} is invalid: {}",
    "W006":
        "No density values provided",
    "W007":
        "No distillation data provided",
    "W008":
        "No reference year provided",
    "W009":
        "Distillation fraction recovered is missing or invalid",
    "W010":
        "Temperature: {} is close to {} -- looks like it could be a K to C conversion error",
}

WARNINGS = {code: (code + ": " + msg) for code, msg in WARNINGS.items()}
