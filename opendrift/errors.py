class WrongMode(Exception):
    def __init__(self, expected_mode, real_mode, msg = None):
        super().__init__(f"Cannot call this function in this mode: {real_mode}, only in: {expected_mode}: {msg}")

class NotCoveredError(Exception):
    pass

class OutsideSpatialCoverageError(NotCoveredError):
    pass

class OutsideTemporalCoverageError(NotCoveredError):
    pass

class VariableNotCoveredError(NotCoveredError):
    pass
