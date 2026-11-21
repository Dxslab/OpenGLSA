from sympy import symbols, Derivative, Function

DERIV_NAME_MAP = {
    'u_fluc':  'uf',
    'v_fluc':  'vf',
    'w_fluc':  'wf',
    'T_fluc':  'Tf',
    'mu_fluc': 'muf',
}

def replace_derivatives(expr):
    replacements = {}

    derivatives = expr.atoms(Derivative)

    for deriv in derivatives:
        func = deriv.expr
        vars_tuple = deriv.variables

        if hasattr(func, 'func'):
            func_name = str(func.func).split('(')[0]
        else:
            func_name = str(func)

        deriv_name = DERIV_NAME_MAP.get(func_name, func_name)

        if len(vars_tuple) == 1:
            var_name = str(vars_tuple[0])
            if var_name == 'xi':
                if func_name in ['Akx', 'Aky', 'Aix', 'Aiy']:
                    symbol_name = f'{func_name}k'
                else:
                    symbol_name = f'd{deriv_name}dk'
            elif var_name == 'eta':
                if func_name in ['Akx', 'Aky', 'Aix', 'Aiy']:
                    symbol_name = f'{func_name}i'
                else:
                    symbol_name = f'd{deriv_name}di'
            else:
                symbol_name = f'd{deriv_name}d{var_name}'
        elif len(vars_tuple) == 2:
            var1, var2 = vars_tuple
            var1_name = str(var1)
            var2_name = str(var2)
            if var1_name == var2_name:
                if var1_name == 'xi':
                    symbol_name = f'd2{deriv_name}dk2'
                elif var1_name == 'eta':
                    symbol_name = f'd2{deriv_name}di2'
                else:
                    symbol_name = f'd2{deriv_name}d{var1_name}2'
            else:
                symbol_name = f'd2{deriv_name}dkdi'
        else:
            var_str = ''.join([f'd{str(vb)}' for vb in vars_tuple])
            symbol_name = f'd{len(vars_tuple)}{deriv_name}{var_str}'

        replacements[deriv] = symbols(symbol_name)

    xi, eta = symbols('xi eta')
    function_names = [
        'Akx', 'Aky', 'Aix', 'Aiy',
        'u', 'v', 'w', 'T', 'mu',
        'u_fluc', 'v_fluc', 'w_fluc', 'T_fluc', 'mu_fluc',
        'x', 'y',
    ]

    for func_name in function_names:
        func = Function(func_name)(xi, eta)
        replacements[func] = symbols(func_name)

    return expr.subs(replacements)