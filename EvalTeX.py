import sympy
import re
from collections import defaultdict

spacing_catching = r'(\\!)|(\\ )|(\\;)|(\\:)|(\\,)|(\\qquad)|(\\quad)|(\\left)|(\\right)'
diff_var_catching = r'((d[a-zA-Z]))'
partial_diff_var_catching = r'(\\partial[a-zA-Z])'
int_var_catching = r'(d[a-zA-Z])'
exponent_catching = r'(\^)'


def unnested(l: 'a nested list') -> [object]:
    if isinstance(l, list):
        if len(l) > 0:
            return unnested(l[0]) + unnested(l[1:])
        else:
            return []
    else:
        return [l]


class ETex:
    '''
    Input:  str (LaTeX expression)
    Output: str (SymPy expression)
    Description:
    The main object which takes one LaTeX expr input and try to parse and translate into a calculable SymPy
    expr. THere are three steps in order to parse a LaTeX expr.
    1. Analyse the total depths/layers within the given LaTeX expr.
            For example. '\int_{1}^{2} \sqrt{\frac{1}{x}} dx' has two depths in total.
                In depth 0, we have ['\int_','^','\sqrt','dx']
                In depth 1, we have ['1','2','\frac']
                In depth 2, we have ['1','x']
    2. Find and catch certain pattern to categorize the raw expr depth by depth.
            in depth 0, the first item '\int_' contains the str '\int' which represents that the input expr is an
            integral. The the underscore of "\int_" show that this integral has a upper and lower bound, therefore it is
            a definite integral.
            Since we know the input expr is a definite integral, then we can analyse this expr in the following way:
                '\int_{lower bound}^{upper bound} (integrand) dx"
            Therefore, we can find the lower bound in depth 1 at the index 0; the upper bound in depth 1 at the index 1;
            And everything except 'dx' in the depth 0 will be counted as the integrand.
    3. Replace each pattern and package them into different parser.
            First we can locate the variable in the term 'dx'. Therefore the variable of this integral is 'x'.
            Then we put everything (upper/lower bound, integrand, variables) in our integral parser, then we finally get:
                'integrate({}, ({}, {}, {}))'.format(integrand, variable, lower, upper)
            Then we can calculate and simplify the raw LaTeX expr by giving SymPy this parsed SymPy expr.
    '''
    @staticmethod
    def convert_frac(num, denum):
        return r'({})/({})'.format(num, denum)

    @staticmethod
    def convert_lim(expr, var, limit):
        return r'limit({}, {}, {})'.format(expr, var, limit)

    @staticmethod
    def convert_int(expr, var, low = None, up = None):
        if low != None and up != None:
            return r'integrate({}, ({}, {}, {}))'.format(expr, var, low, up)    # Definite Integral
        else:
            return r'integrate({}, ({}))'.format(expr, var)                     # Indefinite Integral

    @staticmethod
    def convert_diff(expr, var, power = 1):
        return r'diff({}, {})'.format(expr, ','.join([var] * power))

    @staticmethod
    def convert_partial_diff(expr, var_power_dict):
        starter = r'diff({}, '.format(expr)
        for var, power in var_power_dict.items():
            starter += ','.join([var, repr(power)]) + ','
        starter = starter[:-1] + ')'
        return starter

    def __init__(self, expression: str):
        self.raw_expression = re.sub(spacing_catching, '', expression.replace(" ",""))
        self.depth = defaultdict(list)
        self.max_depth = 0
        self.expr_b4_eval = []
        self.variables = set()
        self.sympy_expr = ''

    def calc(self) -> (set, str):
        self.sympy_expr = self.parse()
        return (self.variables, self.sympy_expr)

    def parse(self) -> str:
        sympy_expr = ''
        self.get_depth()
        if len(self.depth.keys()) > 1:
            self.expr_b4_eval = self.align_depth0_1()
            while self.expr_b4_eval[0]:
                depth0_sec = self.expr_b4_eval[0].pop(0)
                if depth0_sec != '':
                    sympy_expr += self.eval(depth0_sec)
        else:
            sympy_expr = self.eval(self.raw_expression)
        return sympy_expr

    def get_depth(self) -> None:
        has_depth = False
        index_counter = 0
        depth_counter = 0

        for l in self.raw_expression:
            if l == "{":
                has_depth = True
                depth_counter += 1
                self.depth[depth_counter].append(index_counter + 1)
            elif l == "}":
                if depth_counter == 0:
                    raise SyntaxError
                else:
                    self.depth[depth_counter].append(index_counter)
                    depth_counter -= 1
            index_counter += 1

        for key in self.depth.keys():
            if len(self.depth[key]) % 2 != 0:
                raise SyntaxError
            else:
                self.depth[key] = list(zip(self.depth[key][0::2], self.depth[key][1::2]))

        if has_depth:
            self.depth[0].append((0,self.depth[1][0][0] - 1))
            for pair_index in range(0, len(self.depth[1]) - 1):
                if (self.depth[1][pair_index][1] - self.depth[1][pair_index + 1][0] != -2) and \
                   (self.depth[1][pair_index + 1][0] > self.depth[1][pair_index][1]):
                    self.depth[0].append((self.depth[1][pair_index][1] + 1, self.depth[1][pair_index + 1][0] - 1))
            if (self.depth[1][-1][1] + 1) < len(self.raw_expression):
                self.depth[0].append((self.depth[1][-1][1] + 1, len(self.raw_expression)))
            self.max_depth = max(i for i in self.depth.keys())
        else:
            self.depth[0].append((0,len(self.raw_expression)))

    def align_depth0_1(self) -> [list, list]:
        depth_0 = list(filter(lambda x: x != '', unnested([self.raw_expression[pair[0]: pair[1]].split('\\') for pair in self.depth[0]])))
        depth_1 = [self.raw_expression[pair[0]: pair[1]] for pair in self.depth[1]]
        return [depth_0, depth_1]

    def eval(self, depth0_sec) -> str:
        result = ''

        if "frac" in depth0_sec or "dfrac" in depth0_sec:
            var, frac = self.frac_solver(self.expr_b4_eval[1].pop(0), self.expr_b4_eval[1].pop(0))
            self.variables = self.variables.union(var)
            result += frac
        elif "sqrt" in depth0_sec:
            var, int_expr = self.sqrt_solver(depth0_sec)
            self.variables = self.variables.union(var)
            result += int_expr
        elif "int" in depth0_sec:
            var, int_expr = self.int_solver(depth0_sec)
            self.variables = self.variables.union(var)
            result += int_expr
        elif "lim" in depth0_sec:
            var, lim_expr = self.lim_solver(depth0_sec)
            self.variables = self.variables.union(var)
            result += lim_expr
        elif ("sin" in depth0_sec) or ("cos" in depth0_sec) or ("tan" in depth0_sec) \
                or ("sec" in depth0_sec) or ("csc" in depth0_sec) or ("cot" in depth0_sec):
            var, trig_expr = self.trig_func_solver(depth0_sec)
            self.variables = self.variables.union(var)
            result += trig_expr
        elif "^" in depth0_sec:
            var, exponent = self.exponent_solver(depth0_sec)
            self.variables = self.variables.union(var)
            result += exponent
        elif "e" in depth0_sec:
            e_expr = re.sub(r'(e)', 'E', depth0_sec)
            result += e_expr
        elif "pi" in depth0_sec:
            pi_expr = re.sub(r'(\\pi)', '(pi)', depth0_sec)
            result += pi_expr
        elif "infty" in depth0_sec:
            infinity_expr = re.sub(r'(\\infty)', 'oo', depth0_sec)
            result += infinity_expr
        else:
            self.variables = self.variables.union(re.findall(r'([a-zA-Z]+)', depth0_sec))
            result += depth0_sec
        return result

    def frac_solver(self, raw_nom, raw_denom) -> (set, str):
        result = ''
        variables = set()
        if ('d' in raw_nom and 'd' in raw_denom) or ('\\mathrm{d}'in raw_nom and '\\mathrm{d}' in raw_denom):
            return self.derivative_solver(raw_nom, raw_denom)
        elif '\\partial' in raw_nom and '\\partial' in raw_denom:
            return self.partial_derivative_solver(raw_nom, raw_denom)
        if self.max_depth > 1:
            var_in_nominator, nominator = ETex(raw_nom).calc()
            var_in_denominator, denominator = ETex(raw_denom).calc()
            variables = var_in_nominator.union(var_in_denominator)
            result += self.convert_frac(nominator, denominator)
        else:
            result += self.convert_frac(raw_nom, raw_denom)
        return variables, result

    def sqrt_solver(self, raw_sqrt_expr) -> (set, str):
        result = ''
        variables = set()
        if len(raw_sqrt_expr) == 4:
            var, inner_expr = ETex(self.expr_b4_eval[1].pop(0)).calc()
            result += '{}({})'.format(raw_sqrt_expr[0:4], inner_expr)
            variables = var
        return variables, result

    def derivative_solver(self, raw_nom, raw_denom) -> (set, str):
        result = ''
        variables = []
        expr = ''
        p = power = 1
        raw_nom = raw_nom.split('d')
        raw_denom = raw_denom.split('d')

        if len(raw_nom) == 2 and raw_nom[1] == '':
            expr = self.eval(self.expr_b4_eval[0].pop(0))
        elif len(raw_nom) == 2 and raw_nom[1].isalpha():
            expr = raw_nom[1]
        elif len(raw_nom) == 2 and '^{' in raw_nom[1]:
            p = int(re.findall(r'\^{([\d]+)}',raw_nom[1])[0])
            expr = self.eval(self.expr_b4_eval[0].pop(0))
        else:
            raise AssertionError

        if len(raw_denom) == 2:
            if raw_denom[1] == '':
                raise AssertionError
            elif raw_denom[1].isalpha():
                variables.append(raw_denom[1])
            elif '^{' in raw_denom[1]:
                variables.append(raw_denom[1][0: raw_denom[1].find('^{')])
                power = int(re.findall(r'\^{([\d]+)}',raw_denom[1])[0])

        if p == power:
            result += self.convert_diff(expr, variables[0], power)
            return variables, result
        else:
            raise AssertionError


    def partial_derivative_solver(self, raw_nom, raw_denom) -> (set, str):
        result = ''
        variables = []
        expr = ''
        p = 1
        var_power_dict = defaultdict(int)
        raw_nom = raw_nom.split('\\partial')
        raw_denom = raw_denom.split('\\partial')

        if len(raw_nom) == 2 and raw_nom[1] == '':
            expr = self.eval(self.expr_b4_eval[0].pop(0))
        elif len(raw_nom) == 2 and raw_nom[1].isalpha():
            expr = raw_nom[1]
        elif len(raw_nom) == 2 and '^{' in raw_nom[1]:
            p = int(re.findall(r'\^{([\d]+)}',raw_nom[1])[0])
            if raw_nom[1].endswith('}'):
                expr = self.eval(self.expr_b4_eval[0].pop(0))
            else:
                expr = raw_nom[1][raw_nom[1].rfind('}')+1:]
        else:
            raise AssertionError

        if len(raw_denom) >= 2:
            if raw_denom[1] == '':
                raise AssertionError
            for raw_var in raw_denom[1:]:
                if raw_var.isalpha():
                    var_power_dict[raw_var] += 1
                    variables.append(raw_var)
                elif '^{' in raw_var:
                    var = raw_var[0: raw_var.find('^{')]
                    power = int(re.findall(r'\^{([\d]+)}', raw_var)[0])
                    var_power_dict[var] += power

        if p == sum(var_power_dict.values()):
            result += self.convert_partial_diff(expr, var_power_dict)
            return variables, result
        else:
            raise AssertionError


    def lim_solver(self, raw_lim_expr) -> (set, str):
        result =''
        variables = []
        limit = ''

        try:
            raw_var, raw_limit = self.expr_b4_eval[1].pop(0).split('\\to')
            variables += re.findall(r'([a-zA-Z])', raw_var)
            limit += ETex(raw_limit).calc()[1]
        except:
            raise SyntaxError

        lim_expr = self.eval(self.expr_b4_eval[0].pop(0))
        result += self.convert_lim(lim_expr, variables[0], limit)
        return variables, result


    def int_solver(self, raw_int_expr) -> (set, str):
        result = ''
        variables = []
        lower = None
        upper = None

        # Find out the sequence of upper and lower bounds.
        try:
            bound = raw_int_expr[3]
            if bound == "_":
                lower = ETex(self.expr_b4_eval[1].pop(0)).calc()[1]
                upper = ETex(self.expr_b4_eval[1].pop(0)).calc()[1]
                self.expr_b4_eval[0].pop(0)
            elif bound == "^":
                upper = ETex(self.expr_b4_eval[1].pop(0)).calc()[1]
                lower = ETex(self.expr_b4_eval[1].pop(0)).calc()[1]
                self.expr_b4_eval[0].pop(0)
        except IndexError:
            pass

        # Analyse the following expressions.
        try:
            function_expr = self.eval(self.expr_b4_eval[0].pop(0))
        except IndexError:
            function_expr = self.eval(raw_int_expr[4:])
        try:
            # For the situation which 'dx' is at the nominator
            raw_var = re.search(int_var_catching, function_expr).groups()
            var = raw_var[0][1:]
            variables.append(var)
            function_expr = re.sub(int_var_catching, '(1)', function_expr)
        except AttributeError:
            try:
                raw_var = re.search(int_var_catching, self.expr_b4_eval[0][0][0:2]).groups()
                var = raw_var[0][1]
                self.expr_b4_eval[0][0] = self.expr_b4_eval[0][0][2:]
                self.variables.add(var)
            except:
                raise SyntaxError
        result += self.convert_int(function_expr, var, lower, upper)
        return variables, result

    def trig_func_solver(self, raw_trig_expr) -> (set, str):
        result = ''
        variables = set()
        if len(raw_trig_expr) == 3:
            var, inner_expr = ETex(self.expr_b4_eval[1].pop(0)).calc()
            result += '{}({})'.format(raw_trig_expr[0:3], inner_expr)
            variables = variables.union(var)
        elif raw_trig_expr[3] == "^":
            variables, exponent = ETex(self.expr_b4_eval[1].pop(0)).calc()
            var, inner_expr = ETex(self.expr_b4_eval[1].pop(0)).calc()
            result += '{}({})**({})'.format(raw_trig_expr[0:3], inner_expr, exponent)
            variables = variables.union(var)
        return variables, result

    def exponent_solver(self, raw_exponent) -> (set, str):
        result = ''
        variables = set()
        var, exponent = ETex(self.expr_b4_eval[1].pop(0)).calc()
        result += '{}({})'.format(re.sub(exponent_catching, '**', raw_exponent), exponent)
        if "e" in result:
            result = re.sub(r'(e)', 'E', result)
        elif "pi" in result:
            result = re.sub(r'(\\pi)', '(pi)', result)
        variables = variables.union(var)
        for var_in_expr in re.findall(r'([a-zA-Z^]+)', result):
            if var_in_expr not in variables and var_in_expr not in ['E', 'pi']:
                variables.add(var_in_expr)
        return variables, result