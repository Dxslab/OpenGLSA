from sympy import symbols, Function, diff, Rational, simplify, collect, cos, sin, expand
from sympy.printing.fortran import FCodePrinter
import replace_derivatives as rd

xi, eta, z = symbols('xi eta z')
Cp, Pr = symbols('Cp Pr')
wavenumber = symbols('wavenumber')

x   = Function('x')(xi, eta)
y   = Function('y')(xi, eta)
Akx = Function('Akx')(xi, eta)
Aky = Function('Aky')(xi, eta)
Aix = Function('Aix')(xi, eta)
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
    return diff(f, xi) * Akx + diff(f, eta) * Aix

def d_dy(f):
    return diff(f, xi) * Aky + diff(f, eta) * Aiy

def d_dz(f):
    return diff(f, z)


s11b = mu*( 2*d_dx(u)           - 2*Rational(1, 3)*(d_dx(u) + d_dy(v) + v/y + d_dz(w)/y) )
s22b = mu*( 2*d_dy(v)           - 2*Rational(1, 3)*(d_dx(u) + d_dy(v) + v/y + d_dz(w)/y) )
s33b = mu*( 2*d_dz(w)/y + 2*v/y - 2*Rational(1, 3)*(d_dx(u) + d_dy(v) + v/y + d_dz(w)/y) )
s12b = mu*( d_dy(u)   + d_dx(v) )
s13b = mu*( d_dz(u)/y + d_dx(w) )
s23b = mu*( d_dz(v)/y + d_dy(w) - w/y )

s11f = mu*( 2*d_dx(uf)          - 2*Rational(1, 3)*(d_dx(uf)+ d_dy(vf)+ vf/y+ d_dz(wf)/y) ) + \
      muf*( 2*d_dx(u)           - 2*Rational(1, 3)*(d_dx(u) + d_dy(v) + v/y + d_dz(w)/y ) )

s22f = mu*( 2*d_dy(vf)          - 2*Rational(1, 3)*(d_dx(uf)+ d_dy(vf)+ vf/y+ d_dz(wf)/y) ) + \
      muf*( 2*d_dy(v)           - 2*Rational(1, 3)*(d_dx(u) + d_dy(v) + v/y + d_dz(w)/y ) )

s33f = mu*( 2*d_dz(wf)/y+ 2*vf/y- 2*Rational(1, 3)*(d_dx(uf)+ d_dy(vf)+ vf/y+ d_dz(wf)/y) ) + \
      muf*( 2*d_dz(w)/y + 2*v/y - 2*Rational(1, 3)*(d_dx(u) + d_dy(v) + v/y + d_dz(w)/y ) )

s12f = mu*( d_dy(uf)    + d_dx(vf) ) + \
          muf*( d_dy(u) + d_dx(v) )

s13f = mu*( d_dz(uf)/y    + d_dx(wf) ) + \
          muf*( d_dz(u)/y + d_dx(w) )

s23f = mu*( d_dz(vf)/y    + d_dy(wf)    - wf/y ) + \
          muf*( d_dz(v)/y + d_dy(w) - w/y )

b1f = u*s11f + v*s12f + w*s13f + (Cp/Pr)*mu*d_dx(Tf) + \
      uf*s11b    + vf*s12b    + wf*s13b    + (Cp/Pr)*muf*d_dx(T)

b2f = u*s12f + v*s22f + w*s23f + (Cp/Pr)*mu*d_dy(Tf) + \
      uf*s12b    + vf*s22b    + wf*s23b    + (Cp/Pr)*muf*d_dy(T)

b3f = u*s13f + v*s23f + w*s33f + (Cp/Pr)*mu*d_dz(Tf)/y + \
      uf*s13b    + vf*s23b    + wf*s33b    + (Cp/Pr)*muf*d_dz(T)/y


E1 = d_dx(s11f) + d_dy(s12f) + d_dz(s13f)/y + s12f/y
E2 = d_dx(s12f) + d_dy(s22f) + d_dz(s23f)/y + s22f/y - s33f/y
E3 = d_dx(s13f) + d_dy(s23f) + d_dz(s33f)/y + 2*s23f/y
E4 = d_dx(b1f ) + d_dy(b2f ) + d_dz(b3f )/y + b2f/y

metric_dy_dx = diff(y, xi) * Akx + diff(y, eta) * Aix
metric_dy_dy = diff(y, xi) * Aky + diff(y, eta) * Aiy
E1 = E1.subs(metric_dy_dx, 0).subs(metric_dy_dy, 1)
E2 = E2.subs(metric_dy_dx, 0).subs(metric_dy_dy, 1)
E3 = E3.subs(metric_dy_dx, 0).subs(metric_dy_dy, 1)
E4 = E4.subs(metric_dy_dx, 0).subs(metric_dy_dy, 1)

E1 = expand(simplify(E1/cos(wavenumber*z)))
E2 = expand(simplify(E2/cos(wavenumber*z)))
E3 = expand(simplify(E3/sin(wavenumber*z)))
E4 = expand(simplify(E4/cos(wavenumber*z)))

E1_simplified = rd.replace_derivatives(E1)
E2_simplified = rd.replace_derivatives(E2)
E3_simplified = rd.replace_derivatives(E3)
E4_simplified = rd.replace_derivatives(E4)

y_sym = symbols('y')
inv_y = symbols('inv_y')

E1_simplified = E1_simplified.subs(y_sym, 1/inv_y)
E2_simplified = E2_simplified.subs(y_sym, 1/inv_y)
E3_simplified = E3_simplified.subs(y_sym, 1/inv_y)
E4_simplified = E4_simplified.subs(y_sym, 1/inv_y)

def split_y_terms(expr):
    expr_expanded = expand(expr)
    term_reg = expr_expanded.subs(inv_y, 0)
    term_invy = expr_expanded - term_reg
    return term_reg, term_invy

E1_reg, E1_invy = split_y_terms(E1_simplified)
E2_reg, E2_invy = split_y_terms(E2_simplified)
E3_reg, E3_invy = split_y_terms(E3_simplified)
E4_reg, E4_invy = split_y_terms(E4_simplified)

mu_fluc = symbols('mu_fluc')
dmufdk = symbols('dmufdk')
dmufdi = symbols('dmufdi')
mu = symbols('mu')
dmudk = symbols('dmudk')
dmudi = symbols('dmudi')

mu_vars = [mu_fluc, mu, dmufdk, dmudk, dmufdi, dmudi]

E1_reg = collect(E1_reg, mu_vars)
E1_invy = collect(E1_invy, mu_vars)
E2_reg = collect(E2_reg, mu_vars)
E2_invy = collect(E2_invy, mu_vars)
E3_reg = collect(E3_reg, mu_vars)
E3_invy = collect(E3_invy, mu_vars)
E4_reg = collect(E4_reg, mu_vars)
E4_invy = collect(E4_invy, mu_vars)

onethird, twothird, fourthird, eighthird, tenthird = symbols('onethird twothird fourthird eighthird tenthird')
frac_replacements = {
    Rational(10, 3): tenthird,
    Rational(8, 3): eighthird,
    Rational(4, 3): fourthird,
    Rational(2, 3): twothird,
    Rational(1, 3): onethird,
}

E1_reg = E1_reg.subs(frac_replacements)
E1_invy = E1_invy.subs(frac_replacements)
E2_reg = E2_reg.subs(frac_replacements)
E2_invy = E2_invy.subs(frac_replacements)
E3_reg = E3_reg.subs(frac_replacements)
E3_invy = E3_invy.subs(frac_replacements)
E4_reg = E4_reg.subs(frac_replacements)
E4_invy = E4_invy.subs(frac_replacements)

class MyFortranPrinter(FCodePrinter):
    def _wrap_fortran(self, lines):
        return lines

    def _print_Symbol(self, expr):
        name = expr.name
        
        if name in ['xi', 'eta', 'inv_y', 'onethird', 'twothird', 'fourthird', 'eighthird', 'tenthird']:  
            return name
            
        if name in ['Cp', 'Pr', 'wavenumber']: 
            return f"Para%{name}"
            
        return f"{name}(i,j)"
    
    def _print_Integer(self, expr):
        if hasattr(self, '_in_pow') and self._in_pow:
            return super()._print_Integer(expr)
        return f"{expr.p}.0D0"

    def _print_Float(self, expr):
        res = super()._print_Float(expr)
        return res.replace('e', 'D').replace('E', 'D')
        
    def _print_Pow(self, expr):
        self._in_pow = True
        res = super()._print_Pow(expr)
        self._in_pow = False
        return res

    def _print_Add(self, expr):
        terms = expr.as_ordered_terms()
        if not terms:
            return "0.0D0"

        max_line_length = 100
        
        output = self._print(terms[0])
        current_line_length = len(output) + 30
        
        for term in terms[1:]:
            term_str = self._print(term)
            if term_str.startswith('-'):
                separator = " - "
                term_str = term_str[1:].strip()
            else:
                separator = " + "
                
            addition = separator + term_str
            
            if current_line_length + len(addition) > max_line_length:
                output += " " + separator.strip() + " &\n      " + term_str
                current_line_length = 6 + len(term_str)
            else:
                output += addition
                current_line_length += len(addition)
                
        return output

def print_equation(f, expr_reg, expr_invy, eq_name, printer):
    f.write(f"      ! =========================================\n")
    f.write(f"      ! --- {eq_name} Equation Components ---\n")
    f.write(f"      ! =========================================\n")
    
    for part_name, expr in [("reg", expr_reg), ("invy", expr_invy)]:
        f.write(f"      ! --- {part_name} part ---\n")
        
        term_mu_fluc = expr.coeff(mu_fluc)
        term_dmufdk = expr.coeff(dmufdk)
        term_dmufdi = expr.coeff(dmufdi)
        term_dmudi = expr.coeff(dmudi)
        term_dmudk = expr.coeff(dmudk)
        term_mu = expr.coeff(mu)

        f.write(f"      dmufdk_term_{eq_name}_{part_name} = {printer.doprint(term_dmufdk)}\n \n")
        f.write(f"      dmufdi_term_{eq_name}_{part_name} = {printer.doprint(term_dmufdi)}\n \n")
        f.write(f"      muf_term_{eq_name}_{part_name} = {printer.doprint(term_mu_fluc)}\n \n")
        f.write(f"      dmudk_term_{eq_name}_{part_name} = {printer.doprint(term_dmudk)}\n \n")
        f.write(f"      dmudi_term_{eq_name}_{part_name} = {printer.doprint(term_dmudi)}\n \n")
        f.write(f"      mu_term_{eq_name}_{part_name} = {printer.doprint(term_mu)}\n \n")

with open('Cyl_Viscous.inc', 'w', encoding='utf-8') as f:
    printer = MyFortranPrinter({'source_format': 'free', 'standard': 2003})
        
    print_equation(f, E1_reg, E1_invy, "E1", printer)
    print_equation(f, E2_reg, E2_invy, "E2", printer)
    print_equation(f, E3_reg, E3_invy, "E3", printer)
    print_equation(f, E4_reg, E4_invy, "E4", printer)

