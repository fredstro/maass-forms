"""
Base exception hierarchy for the Maass forms ecosystem.

All domain-specific packages (maass-forms-hilbert, maass-forms-klein) should
inherit from these base exceptions to provide a consistent error handling
interface across the ecosystem.

EXAMPLES::

    sage: from maass_form_core.exceptions import MaassFormError, ValidationError
    sage: try:
    ....:     raise ValidationError("bad input", field_name="M", value=-1)
    ....: except MaassFormError as e:
    ....:     print(e)
    bad input
"""

from typing import Any


class MaassFormError(Exception):
    """Base exception for all Maass form packages.

    All package-specific exceptions inherit from this class to provide
    a common exception hierarchy for error handling.
    """

    pass


class InvalidSpectralParameterError(MaassFormError):
    """Raised when a spectral parameter is invalid or out of bounds.

    This exception is raised when:
    - Spectral parameter values are outside expected mathematical bounds
    - Parameter format is incorrect for the mathematical context
    - Complex parameters have invalid real/imaginary parts
    """

    def __init__(self, message: str, parameter_value: Any = None):
        super().__init__(message)
        self.parameter_value = parameter_value


class InvalidSpaceError(MaassFormError):
    """Raised when space configuration is invalid.

    This exception is raised when:
    - MaassFormSpace initialization parameters are invalid
    - Space serialization/deserialization fails
    - Required space properties are missing
    """

    def __init__(self, message: str, space_config: dict | None = None):
        super().__init__(message)
        self.space_config = space_config


class ComputationError(MaassFormError):
    """Raised when mathematical computations fail.

    This exception is raised when:
    - Numerical computations don't converge
    - Mathematical operations encounter singularities
    - Algorithm-specific failures occur
    """

    def __init__(self, message: str, computation_details: dict | None = None):
        super().__init__(message)
        self.computation_details = computation_details


class DatabaseError(MaassFormError):
    """Raised when database operations fail.

    This exception is raised when:
    - MongoDB connection or query failures occur
    - Data serialization/deserialization fails
    - Database validation errors occur
    """

    def __init__(self, message: str, operation: str | None = None):
        super().__init__(message)
        self.operation = operation


class ValidationError(MaassFormError):
    """Raised when input validation fails.

    This exception is raised when:
    - Method parameters don't meet validation criteria
    - Data types are incorrect for mathematical operations
    - Required parameters are missing
    """

    def __init__(self, message: str, field_name: str | None = None, value: Any = None):
        super().__init__(message)
        self.field_name = field_name
        self.value = value
