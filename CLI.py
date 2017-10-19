import EvalTeX
from IPython.display import display_latex, Latex
from sympy import *

funcs = dict()

def ETeX(command, print_latex = True, print_sympy = True):

    expr = EvalTeX.ETex(command).calc()
    if expr[0] != set():
        var(' '.join(expr[0]))
        exec(','.join(expr[0]) + r' = symbols({})'.format(repr(' '.join(expr[0]))))
    try:
        print("Input:=================================================")
        if print_latex:
            print("LaTeX:")
            display_latex(Latex('$$' + command + '$$'))
        if print_sympy:
            print("Sympy:{} \n".format(expr[1]))
        respond = simplify(eval(expr[1], None, funcs))
        print("Output:================================================")
        if print_latex:
            display_latex(Latex('$$' + latex(respond) + '$$'))
        if print_sympy:
            print("Sympy:{} \n".format(respond))
        return respond
    except SyntaxError:
        exec (expr[1], None, funcs)
        return expr[1]
