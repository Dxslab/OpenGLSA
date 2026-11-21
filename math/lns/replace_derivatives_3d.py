from collections import Counter

from sympy import symbols, Derivative, Function

VAR_SUFFIX_MAP = {
    'xi': 'k',
    'eta': 'i',
    'zeta': 's',
}

MIXED_SECOND_SUFFIX_MAP = {
    frozenset(('xi', 'eta')): 'dkdi',
    frozenset(('xi', 'zeta')): 'dkds',
    frozenset(('eta', 'zeta')): 'dids',
}

METRIC_NAMES = {
    'Akx', 'Aky', 'Akz',
    'Aix', 'Aiy', 'Aiz',
    'Asx', 'Asy', 'Asz',
}

DERIV_NAME_MAP = {
    'u_fluc':  'uf',
    'v_fluc':  'vf',
    'w_fluc':  'wf',
    'T_fluc':  'Tf',
    'mu_fluc': 'muf',
}

FUNCTION_NAMES = [
    'Akx', 'Aky', 'Akz',
    'Aix', 'Aiy', 'Aiz',
    'Asx', 'Asy', 'Asz',
    'u', 'v', 'w', 'T', 'mu',
    'u_fluc', 'v_fluc', 'w_fluc', 'T_fluc', 'mu_fluc',
]


def _base_name(func):
    if hasattr(func, 'func'):
        return str(func.func).split('(')[0]
    return str(func)


def _first_derivative_name(func_name, deriv_name, var_name):
    suffix = VAR_SUFFIX_MAP.get(var_name)
    if suffix is None:
        return f'd{deriv_name}d{var_name}'

    if func_name in METRIC_NAMES:
        return f'{func_name}{suffix}'

    return f'd{deriv_name}d{suffix}'


def _second_derivative_name(deriv_name, vars_tuple):
    var_names = [str(var) for var in vars_tuple]
    counts = Counter(var_names)

    if len(counts) == 1:
        var_name = var_names[0]
        suffix = VAR_SUFFIX_MAP.get(var_name, var_name)
        return f'd2{deriv_name}d{suffix}2'

    if len(counts) == 2 and all(count == 1 for count in counts.values()):
        suffix = MIXED_SECOND_SUFFIX_MAP.get(
            frozenset(var_names),
            ''.join(f'd{VAR_SUFFIX_MAP.get(var_name, var_name)}' for var_name in var_names),
        )
        return f'd2{deriv_name}{suffix}'

    var_str = ''.join(f'd{VAR_SUFFIX_MAP.get(var_name, var_name)}' for var_name in var_names)
    return f'd2{deriv_name}{var_str}'


def replace_derivatives(expr):
    replacements = {}

    derivatives = expr.atoms(Derivative)

    for deriv in derivatives:
        func = deriv.expr
        vars_tuple = deriv.variables

        func_name = _base_name(func)

        deriv_name = DERIV_NAME_MAP.get(func_name, func_name)

        if len(vars_tuple) == 1:
            symbol_name = _first_derivative_name(func_name, deriv_name, str(vars_tuple[0]))
        elif len(vars_tuple) == 2:
            symbol_name = _second_derivative_name(deriv_name, vars_tuple)
        else:
            var_str = ''.join(f'd{VAR_SUFFIX_MAP.get(str(var), str(var))}' for var in vars_tuple)
            symbol_name = f'd{len(vars_tuple)}{deriv_name}{var_str}'

        replacements[deriv] = symbols(symbol_name)

    xi, eta, zeta = symbols('xi eta zeta')

    for func_name in FUNCTION_NAMES:
        func = Function(func_name)(xi, eta, zeta)
        replacements[func] = symbols(func_name)

    return expr.subs(replacements)
