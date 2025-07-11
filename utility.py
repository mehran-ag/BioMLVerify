from colorama import Fore, Back, Style, init

import traceback
import sys
import os
import exceptions
import sympy as sp

init( autoreset=True )

color_map = {
    "white": Fore.WHITE,
    "blue": Fore.BLUE,
    "red": Fore.RED,
    "green": Fore.GREEN,
    "yellow": Fore.YELLOW,
    "cyan": Fore.CYAN,
    "magenta": Fore.MAGENTA,
    "light_blue": Fore.LIGHTBLUE_EX,
    "light_cyan": Fore.LIGHTCYAN_EX,
    "light_red": Fore.LIGHTRED_EX,
    "light_green": Fore.LIGHTGREEN_EX,
    "light_yellow": Fore.LIGHTYELLOW_EX,
    "light_magenta": Fore.LIGHTMAGENTA_EX
}

style_map = {
    "normal": Style.NORMAL,
    "bold": Style.BRIGHT,
    "dim": Style.DIM
}

warnings = []


def printer(descrip, text_to_print, descript_color = "white", text_color = "blue", text_style = "normal", u_end="\n"):

    description = f"{color_map.get(descript_color.lower(), Fore.WHITE)}{descrip}"
    text = f"{style_map.get(text_style.lower(), Style.NORMAL)}{color_map.get(text_color.lower(), Fore.BLUE)}{text_to_print}"

    print(f"{description} {text}", end=u_end)

def error_printer(message, error = f"\nERROR: ", error_color = "magenta", error_style="normal", u_end="\n"):

    message = f"{style_map.get(error_style.lower(), Style.NORMAL)}{color_map.get(error_color.lower(), Fore.BLUE)}{message}"

    print(f"{error} {message}", end=u_end)

def message_printer(message, color="white", style = "normal"):

    message = f"{style_map.get(style.lower(), Style.NORMAL)}{color_map.get(color.lower(), Fore.WHITE)}{message}"

    print(f"{message}")

def warning_printer(message, warning = f"\nWARNING: ", color="yellow", style = "normal"):

    message = f"{style_map.get(style.lower(), Style.NORMAL)}{color_map.get(color.lower(), Fore.WHITE)}{message}"

    print(f"{warning}{message}")





def add_warning(message):

    warnings.append(message)

def display_warnings():

    if len(warnings) == 0:
        return
    else:

        print("\n" + "!"*60 + "\n")
        print( "There is/are warning(s) for your compatibility check that need(s) to be investigated for accurate results. The warning(s) is/are listed below:")

        for warning in warnings:

            print( Fore.YELLOW + warning )





def error_handler(e, function=None, print_trace=True):
    
    if function:
        printer("\nAn error has been raised in function: ", function, text_color='magenta')
    else:
        message_printer("\nOperation failed!", color="red")

    if  isinstance(e, sp.SympifyError):
        error_printer("There are variables in the reaction rate equations that are in Python's Keyword list. This prevents Sympify from converting the equation to Sympy equation to be studied!")
        printer("\nFind the expression below:\n" , e)

    elif isinstance(e, exceptions.NotParsable):
        error_printer(e, "\nLIBSBML ERROR: ")

    elif isinstance(e, exceptions.MaxDepth):
        error_printer(e, "\nMAX RECURSION ERROR: ")
        message_printer("There are too many variables in the equation!", color='magenta')


    elif isinstance(e, (TypeError, FileNotFoundError, ValueError, exceptions.NoModel, exceptions.EmptyList, exceptions.NoReverseRateConstant)):
        error_printer(e)

    elif isinstance(e, RuntimeWarning):
        error_printer(e)
    
    elif print_trace:
        tb = traceback.extract_tb(sys.exc_info()[2])
        file_path, lineno, _, _ = tb[-1]
        file_name = os.path.basename(file_path)
        
        error_printer(e, "\nUNEXPECTED ERROR: ")
        printer("Error occurred in file: ", file_name, u_end="")
        printer(", line: ", lineno)
        printer("Error type: ", sys.exc_info()[0].__name__)
        message_printer("Unable to complete the query!", color="red")