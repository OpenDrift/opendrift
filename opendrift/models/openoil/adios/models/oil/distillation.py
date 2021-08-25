'''
    Classes for storing measured values within an Oil record
'''
from dataclasses import dataclass, field

from ..common.utilities import dataclass_to_json, JSON_List
from ..common.measurement import (Time,
                                  Temperature,
                                  MassFraction,
                                  VolumeFraction,
                                  MassOrVolumeFraction,
                                  Concentration,
                                  Unitless,
                                  Dimensionless,
                                  Density,
                                  DynamicViscosity,
                                  KinematicViscosity,
                                  InterfacialTension,
                                  Pressure,
                                  AngularVelocity)

from ..common.validators import EnumValidator
from .validation.warnings import WARNINGS
from .validation.errors import ERRORS


@dataclass_to_json
@dataclass
class DistCut:
    fraction: MassOrVolumeFraction
    vapor_temp: Temperature


class DistCutList(JSON_List):
    item_type = DistCut


@dataclass_to_json
@dataclass
class Distillation:
    type: str = None
    method: str = None
    end_point: Temperature = None
    fraction_recovered: MassOrVolumeFraction = None
    cuts: DistCutList = field(default_factory=DistCutList)

    def validate(self):
        msgs = []
        if self.cuts:  # only need to validate if there are cuts
            msgs.extend(EnumValidator({"mass fraction", "volume fraction"},
                                      ERRORS["E032"],
                                      case_insensitive=True)(self.type))
            if self.fraction_recovered is None:
                msgs.append(WARNINGS["W009"])
            else:
                frac_recov = self.fraction_recovered.converted_to("fraction")
                if frac_recov.value is not None:
                    val = frac_recov.value
                elif frac_recov.max_value is not None:
                    val = frac_recov.max_value

                else:
                    val = None
                    msgs.append(WARNINGS["W009"])
                if val is not None:
                    if not (0.0 <= val <= 1.0):
                        msgs.append(ERRORS["E041"]
                                    .format("distillation fraction recovered",
                                            val))

            for cut in self.cuts:
                frac = cut.fraction.converted_to('fraction').value
                if not (0.0 <= frac <= 1.0):
                    msgs.append(ERRORS["E041"]
                                .format("distillation fraction", frac))
                vt = cut.vapor_temp.convert_to('C').value
                if vt < -100.0:
                    t = f"{cut.vapor_temp.value:.2f} {cut.vapor_temp.unit}"
                    msgs.append(ERRORS["E040"]
                                .format("distillation vapor temp", t))

        return msgs
