ERRORS = {
    # E01* : meta data related
    "E010": "Record has no oil_id: every record must have an ID",
    "E011": ("Record has invalid oil_id: every record must have a valid ID. "
             "{} is not valid"),
    "E012": "Reference year: {} is not a valid year (between {} and {})",
    "E013": "Review Status: {} is not valid, must be one of {}",
    # E03* -- physical properties related
    "E030": "Oils must have an API",
    "E031": "No Properties data at all",
    "E032": 'Distillation type is "{}", it must be one of: {}',

    # E04* -- units or values
    "E040": "Value for {}: {} is out of range: unit error?",
    "E041": "Value for {}: {} must be between 0 and 1",
    "E042": "Must have a value for {}",
    "E043": "API, {} does not match density at 60F. API should be: {:.1f}",

    # E05* -- duplicates, etc
    "E050": "Duplicate {} in {}",
}

ERRORS = {code: (code + ": " + msg) for code, msg in ERRORS.items()}
