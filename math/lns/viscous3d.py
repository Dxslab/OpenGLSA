from pathlib import Path

from sympy import symbols, Function, diff, Rational, simplify, collect, expand
from sympy.printing.fortran import FCodePrinter
import replace_derivatives_3d as rd

xi, eta, zeta = symbols('xi eta zeta')
Cp, Pr = symbols('Cp Pr')

Akx = Function('Akx')(xi, eta, zeta)
Aky = Function('Aky')(xi, eta, zeta)
Akz = Function('Akz')(xi, eta, zeta)
Aix = Function('Aix')(xi, eta, zeta)
Aiy = Function('Aiy')(xi, eta, zeta)
Aiz = Function('Aiz')(xi, eta, zeta)
Asx = Function('Asx')(xi, eta, zeta)
Asy = Function('Asy')(xi, eta, zeta)
Asz = Function('Asz')(xi, eta, zeta)

u = Function('u')(xi, eta, zeta)
v = Function('v')(xi, eta, zeta)
w = Function('w')(xi, eta, zeta)
T = Function('T')(xi, eta, zeta)
mu = Function('mu')(xi, eta, zeta)

u_fluc = Function('u_fluc')(xi, eta, zeta)
v_fluc = Function('v_fluc')(xi, eta, zeta)
w_fluc = Function('w_fluc')(xi, eta, zeta)
T_fluc = Function('T_fluc')(xi, eta, zeta)
mu_fluc = Function('mu_fluc')(xi, eta, zeta)

def d_dx(f):
    return diff(f, xi) * Akx + diff(f, eta) * Aix + diff(f, zeta) * Asx

def d_dy(f):
    return diff(f, xi) * Aky + diff(f, eta) * Aiy + diff(f, zeta) * Asy

def d_dz(f):
    return diff(f, xi) * Akz + diff(f, eta) * Aiz + diff(f, zeta) * Asz

# Baseflow
s11b = mu*( 2*d_dx(u) - 2*Rational(1, 3)*(d_dx(u) + d_dy(v) + d_dz(w)) )
s22b = mu*( 2*d_dy(v) - 2*Rational(1, 3)*(d_dx(u) + d_dy(v) + d_dz(w)) )
s33b = mu*( 2*d_dz(w) - 2*Rational(1, 3)*(d_dx(u) + d_dy(v) + d_dz(w)) )
s12b = mu*( d_dy(u) + d_dx(v) )
s13b = mu*( d_dz(u) + d_dx(w) )
s23b = mu*( d_dz(v) + d_dy(w) )

# flucflow
s11f = mu*( 2*d_dx(u_fluc) - 2*Rational(1, 3)*(d_dx(u_fluc) + d_dy(v_fluc) + d_dz(w_fluc)) ) + mu_fluc*( 2*d_dx(u) - 2*Rational(1, 3)*(d_dx(u) + d_dy(v) + d_dz(w)) )
s22f = mu*( 2*d_dy(v_fluc) - 2*Rational(1, 3)*(d_dx(u_fluc) + d_dy(v_fluc) + d_dz(w_fluc)) ) + mu_fluc*( 2*d_dy(v) - 2*Rational(1, 3)*(d_dx(u) + d_dy(v) + d_dz(w)) )
s33f = mu*( 2*d_dz(w_fluc) - 2*Rational(1, 3)*(d_dx(u_fluc) + d_dy(v_fluc) + d_dz(w_fluc)) ) + mu_fluc*( 2*d_dz(w) - 2*Rational(1, 3)*(d_dx(u) + d_dy(v) + d_dz(w)) )
s12f = mu*( d_dy(u_fluc) + d_dx(v_fluc) ) + mu_fluc*( d_dy(u) + d_dx(v) )
s13f = mu*( d_dz(u_fluc) + d_dx(w_fluc) ) + mu_fluc*( d_dz(u) + d_dx(w) )
s23f = mu*( d_dz(v_fluc) + d_dy(w_fluc) ) + mu_fluc*( d_dz(v) + d_dy(w) )

b1f = u*s11f + v*s12f + w*s13f + (Cp/Pr)*mu_fluc*d_dx(T) + u_fluc*s11b  + v_fluc*s12b + w_fluc*s13b + (Cp/Pr)*mu*d_dx(T_fluc)
b2f = u*s12f + v*s22f + w*s23f + (Cp/Pr)*mu_fluc*d_dy(T) + u_fluc*s12b  + v_fluc*s22b + w_fluc*s23b + (Cp/Pr)*mu*d_dy(T_fluc)
b3f = u*s13f + v*s23f + w*s33f + (Cp/Pr)*mu_fluc*d_dz(T) + u_fluc*s13b  + v_fluc*s23b + w_fluc*s33b + (Cp/Pr)*mu*d_dz(T_fluc)

E1 = d_dx(s11f) + d_dy(s12f) + d_dz(s13f)
E2 = d_dx(s12f) + d_dy(s22f) + d_dz(s23f)
E3 = d_dx(s13f) + d_dy(s23f) + d_dz(s33f)
E4 = d_dx(b1f ) + d_dy(b2f ) + d_dz(b3f )

E1 = expand(simplify(E1))
E2 = expand(simplify(E2))
E3 = expand(simplify(E3))
E4 = expand(simplify(E4))

E1_simplified = rd.replace_derivatives(E1)
E2_simplified = rd.replace_derivatives(E2)
E3_simplified = rd.replace_derivatives(E3)
E4_simplified = rd.replace_derivatives(E4)

mu_fluc_sym = symbols('mu_fluc')
dmufdk = symbols('dmufdk')
dmufdi = symbols('dmufdi')
dmufds = symbols('dmufds')
mu_sym = symbols('mu')
dmudk = symbols('dmudk')
dmudi = symbols('dmudi')
dmuds = symbols('dmuds')

mu_vars = [mu_fluc_sym, mu_sym, dmufdk, dmudk, dmufdi, dmudi, dmufds, dmuds]

E1_simplified = collect(E1_simplified, mu_vars)
E2_simplified = collect(E2_simplified, mu_vars)
E3_simplified = collect(E3_simplified, mu_vars)
E4_simplified = collect(E4_simplified, mu_vars)

onethird, twothird, fourthird, eighthird, tenthird = symbols('onethird twothird fourthird eighthird tenthird')
frac_replacements = {
    Rational(10, 3): tenthird,
    Rational(8, 3): eighthird,
    Rational(4, 3): fourthird,
    Rational(2, 3): twothird,
    Rational(1, 3): onethird,
}

E1_simplified = E1_simplified.subs(frac_replacements)
E2_simplified = E2_simplified.subs(frac_replacements)
E3_simplified = E3_simplified.subs(frac_replacements)
E4_simplified = E4_simplified.subs(frac_replacements)

class MyFortranPrinter(FCodePrinter):
    def _wrap_fortran(self, lines):
        return lines

    def _print_Symbol(self, expr):
        name = expr.name
        
        if name in ['xi', 'eta', 'zeta', 'onethird', 'twothird', 'fourthird', 'eighthird', 'tenthird']:  
            return name
            
        if name in ['Cp', 'Pr']: 
            return f"Para%{name}"
            
        return f"{name}(i,j,k)"
    
    def _print_Integer(self, expr):
        # 确保幂次方的数字保持整数
        if hasattr(self, '_in_pow') and self._in_pow:
            return super()._print_Integer(expr)
        return f"{expr.p}.0D0"

    def _print_Float(self, expr):
        # 处理浮点数
        res = super()._print_Float(expr)
        return res.replace('e', 'D').replace('E', 'D')
        
    def _print_Pow(self, expr):
        # 给求幂运算加个标记，保证里面的指数不会变成 2.0D0
        self._in_pow = True
        res = super()._print_Pow(expr)
        self._in_pow = False
        return res

    def _print_Add(self, expr):
        terms = expr.as_ordered_terms()
        if not terms:
            return "0.0D0"

        max_line_length = 100  # 设置你想要的每一行的最大允许长度
        
        # 初始第一项
        output = self._print(terms[0])
        current_line_length = len(output) + 20 # 预估前缀类似 "      mu_term1 = " 的长度
        
        # 处理后续项
        for term in terms[1:]:
            term_str = self._print(term)
            
            # 决定符号和项的字符串
            if term_str.startswith('-'):
                separator = " - "
                term_str = term_str[1:].strip() # 移除负号，避免重复
            else:
                separator = " + "
                
            addition = separator + term_str
            
            # 判断加上这项后是否超长
            if current_line_length + len(addition) > max_line_length:
                # 换行
                output += " " + separator.strip() + " &\n      " + term_str
                current_line_length = 6 + len(term_str) # "      & " 的长度加上项长度
            else:
                # 不换行，继续连在后面
                output += addition
                current_line_length += len(addition)
                
        return output
    
def print_equation(f, expr, n, printer):
    term_mu_fluc = expr.coeff(mu_fluc_sym)
    term_dmufdk = expr.coeff(dmufdk)
    term_dmufdi = expr.coeff(dmufdi)
    term_dmufds = expr.coeff(dmufds)
    term_dmudi = expr.coeff(dmudi)
    term_dmudk = expr.coeff(dmudk)
    term_dmuds = expr.coeff(dmuds)
    term_mu = expr.coeff(mu_sym)

    f.write(f"      dmufdk_term({n}) = {printer.doprint(term_dmufdk)}\n \n")
    f.write(f"      dmufdi_term({n}) = {printer.doprint(term_dmufdi)}\n \n")
    f.write(f"      dmufds_term({n}) = {printer.doprint(term_dmufds)}\n \n")
    f.write(f"      muf_term({n}) = {printer.doprint(term_mu_fluc)}\n \n")
    f.write(f"      dmudk_term({n}) = {printer.doprint(term_dmudk)}\n \n")
    f.write(f"      dmudi_term({n}) = {printer.doprint(term_dmudi)}\n \n")
    f.write(f"      dmuds_term({n}) = {printer.doprint(term_dmuds)}\n \n")
    f.write(f"      mu_term({n}) = {printer.doprint(term_mu)}\n \n")


output_path = Path(__file__).with_name('Viscous3d.inc')

with open(output_path, 'w', encoding='utf-8') as f:
    printer = MyFortranPrinter({'source_format': 'free', 'standard': 2003})
        
    print_equation(f, E1_simplified, 1, printer)
    print_equation(f, E2_simplified, 2, printer)
    print_equation(f, E3_simplified, 3, printer)
    print_equation(f, E4_simplified, 4, printer)
