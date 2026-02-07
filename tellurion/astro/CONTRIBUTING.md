# Contributing to AstroPy: Quantity Utilities

This document outlines the improvements made to prepare `quantity_utils.py` for potential contribution to AstroPy.

## Changes Made

### 1. Code Style & Naming Conventions

**Before:**
- `isscalar` → `is_scalar`
- `unitlookup` → `unit_lookup`
- Inconsistent function naming

**After:**
- All parameters use snake_case
- Function names are descriptive and follow PEP 8
- Consistent naming patterns throughout

### 2. Type Hints

All functions now have complete type hints:
```python
def make_quantity(
    value: Union[float, np.ndarray, u.Quantity, Dict[str, Any]],
    unit: Optional[Union[u.Unit, str, Dict[str, Union[u.Unit, str]]]] = None,
    is_scalar: bool = False,
    unit_lookup: Optional[Dict[str, Union[u.Unit, str]]] = None
) -> u.Quantity:
```

### 3. Documentation

**NumPy-style docstrings** for all functions with:
- One-line summary
- Extended description
- Parameters section with types
- Returns section
- Raises section for exceptions
- Examples section with runnable code
- Notes section for important details

### 4. Function Renaming

For better clarity and AstroPy consistency:
- `hstack()` → `hstack_quantities()`
- `vstack()` → `vstack_quantities()`
- `changeunits()` → `change_units()`
- `quantity.to_dict()` → `quantity_to_dict(quantity)` (removed monkey-patching)
- `quantity.to_array()` → `quantity_to_array(quantity)` (removed monkey-patching)
- `quantity.tovector()` → removed (converted to internal helper if needed)

### 5. Removed Monkey-Patching

**Before:**
```python
u.Quantity.to_dict = lambda self: {...}
u.Quantity.to_array = lambda self: {...}
u.Quantity.tovector = lambda self: {...}
```

**After:**
These are now standalone functions that take a Quantity as an argument. This is better practice and more maintainable.

### 6. Error Handling

Improved error messages with context:
```python
raise ValueError(
    f"Keys in value dict {set(value.keys())} must match "
    f"keys in unit dict {set(unit.keys())}"
)
```

### 7. Testing

Comprehensive test suite with:
- Unit tests for each function
- Edge case testing
- Error condition testing
- Integration tests
- Uses pytest framework (AstroPy standard)
- Organized into test classes by function

## Test Coverage

Run tests with:
```bash
pytest test_quantity_utils.py -v
```

Test coverage includes:
- ✅ All public functions
- ✅ Scalar and array inputs
- ✅ Regular and structured quantities
- ✅ Unit conversions
- ✅ Stacking operations
- ✅ Error conditions
- ✅ Edge cases
- ✅ Integration workflows

## Next Steps for AstroPy Contribution

1. **Open an Issue**: Start a discussion on AstroPy's GitHub about adding structured quantity utilities
2. **Write an APE**: Consider drafting an AstroPy Proposal for Enhancement if this requires API changes
3. **Integration**: Work with AstroPy maintainers to integrate into `astropy.units`
4. **Performance**: Profile and optimize if needed
5. **Documentation**: Add to AstroPy's documentation
6. **Continuous Integration**: Ensure tests run in AstroPy's CI pipeline

## API Changes Summary

If you're migrating from the old `quant.py`:

| Old | New |
|-----|-----|
| `make_quantity(val, unit, isscalar=True)` | `make_quantity(val, unit, is_scalar=True)` |
| `changeunits(q, unitlookup)` | `change_units(q, unit_lookup)` |
| `_make_structured_quantity(q, name, scalar)` | `make_structured_quantity(q, name, is_scalar)` |
| `hstack(sqs)` | `hstack_quantities(quantities)` |
| `vstack(sqs)` | `vstack_quantities(quantities)` |
| `q.to_dict()` | `quantity_to_dict(q)` |
| `q.to_array()` | `quantity_to_array(q)` |
| `q.tovector()` | Removed (use internal helper if needed) |

## Module Structure

```
quantity_utils.py
├── Build quantities
│   ├── make_quantity()
│   └── change_units()
├── Build structured quantities
│   ├── make_structured_quantity()
│   ├── hstack_quantities()
│   └── vstack_quantities()
└── Convert structured and unstructured quantities
    ├── quantity_to_dict()
    └── quantity_to_array()
```

## License

When contributing to AstroPy, this code will be licensed under the BSD 3-Clause License to match AstroPy's license.