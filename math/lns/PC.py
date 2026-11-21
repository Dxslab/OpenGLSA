from sympy import symbols, Function, diff, Rational, simplify, collect, cos, sin, expand
from sympy.printing.fortran import FCodePrinter
import replace_derivatives as rd

xi, eta, z = symbols('xi eta z')
wavenumber = symbols('wavenumber')

x   = Function('x')(xi, eta)
y   = Function('y')(xi, eta)
Akx = Function('Akx')(xi, eta)
Aiy = Function('Aiy')(xi, eta)

u = Function('u')(xi, eta)
v = Function('v')(xi, eta)
w = 0
T = Function('T')(xi, eta)
mu = Function('mu')(xi, eta)

u_fluc = Function('u_fluc')(xi, eta)
v_fluc = Function('v_fluc')(xi, eta)
w_fluc = Function('w_fluc')(xi, eta)
T_fluc = Function('T_fluc')(xi, eta)
mu_fluc = Function('mu_fluc')(xi, eta)

uf = u_fluc * cos(wavenumber * z)
vf = v_fluc * cos(wavenumber * z)
wf = w_fluc * sin(wavenumber * z)
Tf = T_fluc * cos(wavenumber * z)
muf = mu_fluc * cos(wavenumber * z)

def d_dx(f):
    return diff(f, xi) * Akx
def d_dy(f):
    return diff(f, eta) * Aiy

def d_dz(f):
    return diff(f, z)

# Baseflow
s11b = mu*( 2*d_dx(u) - 2*Rational(1, 3)*(d_dx(u) + d_dy(v)) )
s22b = mu*( 2*d_dy(v) - 2*Rational(1, 3)*(d_dx(u) + d_dy(v)) )
s33b = mu*( 2*d_dz(w) - 2*Rational(1, 3)*(d_dx(u) + d_dy(v)) )
s12b = mu*( d_dy(u) + d_dx(v) )
s13b = mu*( d_dz(u) + d_dx(w) )
s23b = mu*( d_dz(v) + d_dy(w) )

# flucflow
s11f = mu*( 2*d_dx(uf) - 2*Rational(1, 3)*(d_dx(uf) + d_dy(vf) + d_dz(wf)) ) + muf*( 2*d_dx(u) - 2*Rational(1, 3)*(d_dx(u) + d_dy(v) + d_dz(w)) )
s22f = mu*( 2*d_dy(vf) - 2*Rational(1, 3)*(d_dx(uf) + d_dy(vf) + d_dz(wf)) ) + muf*( 2*d_dy(v) - 2*Rational(1, 3)*(d_dx(u) + d_dy(v) + d_dz(w)) )
s33f = mu*( 2*d_dz(wf) - 2*Rational(1, 3)*(d_dx(uf) + d_dy(vf) + d_dz(wf)) ) + muf*( 2*d_dz(w) - 2*Rational(1, 3)*(d_dx(u) + d_dy(v) + d_dz(w)) )
s12f = mu*( d_dy(uf) + d_dx(vf) ) + muf*( d_dy(u) + d_dx(v) )
s13f = mu*( d_dz(uf) + d_dx(wf) ) + muf*( d_dz(u) + d_dx(w) )
s23f = mu*( d_dz(vf) + d_dy(wf) ) + muf*( d_dz(v) + d_dy(w) )


E1 = d_dx(s11f) + d_dy(s12f) + d_dz(s13f)
E2 = d_dx(s12f) + d_dy(s22f) + d_dz(s23f)

E1 = expand(simplify(E1/cos(wavenumber*z))/Akx)
E2 = expand(simplify(E2/cos(wavenumber*z))/Aiy)

E1_simplified = rd.replace_derivatives(E1)
E2_simplified = rd.replace_derivatives(E2)

mu_fluc = symbols('mu_fluc')
dmufdk = symbols('dmufdk')
dmufdi = symbols('dmufdi')
mu = symbols('mu')
dmudk = symbols('dmudk')
dmudi = symbols('dmudi')

mu_vars = [mu_fluc, mu, dmufdk, dmudk, dmufdi, dmudi]

E1_simplified = collect(E1_simplified, mu_vars)
E2_simplified = collect(E2_simplified, mu_vars)

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

class MyFortranPrinter(FCodePrinter):
    def _wrap_fortran(self, lines):
        return lines

    def _print_Symbol(self, expr):
        name = expr.name
        
        if name in ['xi', 'eta', 'onethird', 'twothird', 'fourthird', 'eighthird', 'tenthird']:  
            return name
            
        if name in ['Cp', 'Pr', 'wavenumber']: 
            return f"Para%{name}"
            
        return f"{name}(i,j)"
    
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
    term_mu_fluc = expr.coeff(mu_fluc)
    term_dmufdk = expr.coeff(dmufdk)
    term_dmufdi = expr.coeff(dmufdi)
    term_dmudi = expr.coeff(dmudi)
    term_dmudk = expr.coeff(dmudk)
    term_mu = expr.coeff(mu)

    f.write(f"      dmufdk_term{n} = {printer.doprint(term_dmufdk)}\n \n")
    f.write(f"      dmufdi_term{n} = {printer.doprint(term_dmufdi)}\n \n")
    f.write(f"      mu_fluc_term{n} = {printer.doprint(term_mu_fluc)}\n \n")
    f.write(f"      dmudk_term{n} = {printer.doprint(term_dmudk)}\n \n")
    f.write(f"      dmudi_term{n} = {printer.doprint(term_dmudi)}\n \n")
    f.write(f"      mu_term{n} = {printer.doprint(term_mu)}\n \n")

with open('bc_pc.inc', 'w', encoding='utf-8') as f:
    printer = MyFortranPrinter({'source_format': 'free', 'standard': 2003})
        
    print_equation(f, E1_simplified, 1, printer)
    print_equation(f, E2_simplified, 2, printer)
