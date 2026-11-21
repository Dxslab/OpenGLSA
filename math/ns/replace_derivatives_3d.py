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
        
        if len(vars_tuple) == 1:
            var_name = str(vars_tuple[0])
            if var_name == 'xi':
                if func_name in ['Akx', 'Aky', 'Akz', 'Aix', 'Aiy', 'Aiz', 'Asx', 'Asy', 'Asz']:
                    symbol_name = f'{func_name}k'
                else:
                    symbol_name = f'd{func_name}dk'
            elif var_name == 'eta':
                if func_name in ['Akx', 'Aky', 'Akz', 'Aix', 'Aiy', 'Aiz', 'Asx', 'Asy', 'Asz']:
                    symbol_name = f'{func_name}i'
                else:
                    symbol_name = f'd{func_name}di'
            elif var_name == 'zeta':
                if func_name in ['Akx', 'Aky', 'Akz', 'Aix', 'Aiy', 'Aiz', 'Asx', 'Asy', 'Asz']:
                    symbol_name = f'{func_name}s'
                else:
                    symbol_name = f'd{func_name}ds'
            else:
                symbol_name = f'd{func_name}d{var_name}'
        elif len(vars_tuple) == 2:
            var1, var2 = vars_tuple
            var1_name = str(var1)
            var2_name = str(var2)
            if var1_name == var2_name:
                if var1_name == 'xi':
                    symbol_name = f'd2{func_name}dk2'
                elif var1_name == 'eta':
                    symbol_name = f'd2{func_name}di2'
                elif var1_name == 'zeta':
                    symbol_name = f'd2{func_name}ds2'
                else:
                    symbol_name = f'd2{func_name}d{var1_name}2'
            else:
                suffix = MIXED_SECOND_SUFFIX_MAP.get(
                    frozenset((var1_name, var2_name)),
                    f'd{VAR_SUFFIX_MAP.get(var1_name, var1_name)}d{VAR_SUFFIX_MAP.get(var2_name, var2_name)}'
                )
                symbol_name = f'd2{func_name}{suffix}'
        else:
            var_str = ''.join([f'd{str(vb)}' for vb in vars_tuple])
            symbol_name = f'd{len(vars_tuple)}{func_name}{var_str}'
        
        replacements[deriv] = symbols(symbol_name)
    
    xi, eta, zeta = symbols('xi eta zeta')
    function_names = ['Akx', 'Aky', 'Akz', 'Aix', 'Aiy', 'Aiz', 'Asx', 'Asy', 'Asz', 'u', 'v', 'w', 'T', 'mu']
    
    for func_name in function_names:
        func = Function(func_name)(xi, eta, zeta)
        replacements[func] = symbols(func_name)
    
    return expr.subs(replacements)
