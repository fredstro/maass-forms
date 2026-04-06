"""
Custom exceptions for the maass-forms-klein package.

This module defines package-specific exception classes for better error handling
and debugging throughout the maass-forms-klein mathematical framework.
"""

from typing import Any, Optional


class KnotMaassError(Exception):
    """Base exception for maass-forms-klein package.
    
    All package-specific exceptions inherit from this class to provide
    a common exception hierarchy for error handling.
    """
    pass


class InvalidSpectralParameterError(KnotMaassError):
    """Raised when a spectral parameter is invalid or out of bounds.
    
    This exception is raised when:
    - Spectral parameter values are outside expected mathematical bounds
    - Parameter format is incorrect for the mathematical context
    - Complex parameters have invalid real/imaginary parts
    """

    def __init__(self, message: str, parameter_value: Any = None):
        super().__init__(message)
        self.parameter_value = parameter_value


class InvalidSpaceError(KnotMaassError):
    """Raised when space configuration is invalid.
    
    This exception is raised when:
    - KleinianMaassFormSpace initialization parameters are invalid
    - Space serialization/deserialization fails
    - Required space properties are missing
    """

    def __init__(self, message: str, space_config: Optional[dict] = None):
        super().__init__(message)
        self.space_config = space_config


class InvalidGroupError(KnotMaassError):
    """Raised when Kleinian group configuration is invalid.
    
    This exception is raised when:
    - Group generators are not valid matrices
    - Group construction fails mathematical validation
    - Group parameters are inconsistent
    """

    def __init__(self, message: str, group_data: Any = None):
        super().__init__(message)
        self.group_data = group_data


class ComputationError(KnotMaassError):
    """Raised when mathematical computations fail.
    
    This exception is raised when:
    - Numerical computations don't converge
    - Mathematical operations encounter singularities
    - Algorithm-specific failures occur
    """

    def __init__(self, message: str, computation_details: Optional[dict] = None):
        super().__init__(message)
        self.computation_details = computation_details


class DatabaseError(KnotMaassError):
    """Raised when database operations fail.
    
    This exception is raised when:
    - MongoDB connection or query failures occur
    - Data serialization/deserialization fails
    - Database validation errors occur
    """

    def __init__(self, message: str, operation: Optional[str] = None):
        super().__init__(message)
        self.operation = operation


class ValidationError(KnotMaassError):
    """Raised when input validation fails.
    
    This exception is raised when:
    - Method parameters don't meet validation criteria
    - Data types are incorrect for mathematical operations
    - Required parameters are missing
    """

    def __init__(self, message: str, field_name: Optional[str] = None, value: Any = None):
        super().__init__(message)
        self.field_name = field_name
        self.value = value
