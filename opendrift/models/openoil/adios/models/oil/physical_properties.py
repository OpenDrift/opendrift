"""
Main class that represents an oil record.

This maps to the JSON used in the DB

Having a Python class makes it easier to write importing, validating etc, code.
"""
from dataclasses import dataclass, field

from .validation.errors import ERRORS

from ..common.utilities import dataclass_to_json, JSON_List

from ..common.measurement import (Temperature,
                                  Density,
                                  DynamicViscosity,
                                  KinematicViscosity,
                                  AngularVelocity,
                                  InterfacialTension)


class RefTempList:
    """
    mixin for all classes that are a list of points with
    reference temperatures
    """

    def validate(self):
        """
        validater for anything that has a list of reference temps

        e.g. density and viscosity

        For viscosity it checks for shear rate as well.
        """
        points_list = self
        data_str = self.__class__.__name__
        msgs = []
        # check for odd temperatures
        for pt in points_list:
            if pt.ref_temp is None:
                msgs.append(ERRORS["E042"]
                            .format(data_str + " reference temp"))
                return msgs

            temp = pt.ref_temp.converted_to('C').value
            if temp is None:
                msgs.append(ERRORS["E042"]
                            .format(data_str + " reference temp"))
                return msgs
            if temp < -100.0:  # arbitrary, but should catch K/C confusion
                t = f"{pt.ref_temp.value:.2f} {pt.ref_temp.unit}"
                msgs.append(ERRORS["E040"].format(data_str, t))

        # check for duplicate temp/shear_rate combos
        temps = []
        for p in points_list:
            temp = p.ref_temp.converted_to('K').value
            try:
                temp = temp + p.shear_rate.value
            except (TypeError, AttributeError):
                pass
            temps.append(temp)
        temps.sort()
        diff = (abs(t2 - t1) for t1, t2 in zip(temps[1:], temps[:1]))
        for d in diff:
            if d < 1e-3:
                msgs.append(ERRORS["E050"].format("Temperatures", data_str))

        return msgs


@dataclass_to_json
@dataclass
class DensityPoint:
    density: Density = None
    ref_temp: Temperature = None
    method: str = None


class DensityList(JSON_List, RefTempList):
    item_type = DensityPoint


@dataclass_to_json
@dataclass
class DynamicViscosityPoint:
    viscosity: DynamicViscosity = None
    ref_temp: Temperature = None
    shear_rate: AngularVelocity = None
    method: str = None


class DynamicViscosityList(JSON_List, RefTempList):
    item_type = DynamicViscosityPoint


@dataclass_to_json
@dataclass
class KinematicViscosityPoint:
    viscosity: KinematicViscosity = None
    ref_temp: Temperature = None
    shear_rate: AngularVelocity = None
    method: str = None


class KinematicViscosityList(JSON_List, RefTempList):
    item_type = KinematicViscosityPoint


@dataclass_to_json
@dataclass
class PourPoint:
    measurement: Temperature = None
    method: str = None


@dataclass_to_json
@dataclass
class FlashPoint:
    measurement: Temperature = None
    method: str = None


@dataclass_to_json
@dataclass
class InterfacialTensionPoint:
    tension: InterfacialTension
    ref_temp: Temperature
    # interface: str = None  # the interface is given in the attribute name
    method: str = None


class InterfacialTensionList(JSON_List, RefTempList):
    item_type = InterfacialTensionPoint


@dataclass_to_json
@dataclass
class PhysicalProperties:
    pour_point: PourPoint = None
    flash_point: FlashPoint = None

    densities: DensityList = field(default_factory=DensityList)
    kinematic_viscosities: KinematicViscosityList = field(default_factory=KinematicViscosityList)
    dynamic_viscosities: DynamicViscosityList = field(default_factory=DynamicViscosityList)

    interfacial_tension_air: InterfacialTensionList = field(default_factory=InterfacialTensionList)
    interfacial_tension_water: InterfacialTensionList = field(default_factory=InterfacialTensionList)
    interfacial_tension_seawater: InterfacialTensionList = field(default_factory=InterfacialTensionList)
