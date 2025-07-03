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


def printer(descrip, text_to_print, descript_color = "white", text_color = "blue", text_style = "normal"):

    description = f"{color_map.get(descript_color.lower(), Fore.WHITE)}{descrip}"
    text = f"{style_map.get(text_style.lower(), Style.NORMAL)}{color_map.get(text_color.lower(), Fore.BLUE)}{text_to_print}"

    print(f"{description} {text}")

def error_printer(message, error, message_color = 'white', error_color = "red", error_style="normal", u_end="\n"):

    message = f"{color_map.get(message_color.lower(), Fore.WHITE)}{message}"
    error = f"{style_map.get(error_style.lower(), Style.NORMAL)}{color_map.get(error_color.lower(), Fore.BLUE)}{error}"

    print(f"{message} {error}", end=u_end)

def message_printer(message, color="yellow", style = "bold"):

    message = f"{style_map.get(style.lower(), Style.NORMAL)}{color_map.get(color.lower(), Fore.WHITE)}{message}"

    print(f"{message}")





def add_warning(message):

    warnings.append(message)

def display_warnings():

    if len(warnings) == 0:
        return
    else:

        print("\n" + "!"*60)
        print( "There is/are warning(s) for your compatibility check that need(s) to be investigated for accurate results. The warning(s) is/are listed below:")

        for warning in warnings:

            print( Fore.LIGHTCYAN_EX + warning )





def error_handler(e, function=None, print_trace=True):
    
    if function:
        printer("\nAn error has been raised in function: ", function, text_color='white')
    else:
        message_printer("\nOperation failed!", color="red")

    if isinstance(e, (TypeError, FileNotFoundError, ValueError, exceptions.NoModel, exceptions.EmptyList, exceptions.NoReverseRateConstant)):
        printer("ERROR: ", e, text_color='cyan')

    elif isinstance(e, sp.SympifyError):
        error_printer("Sympify Error: ", e)
        message_printer("Equation couldn't be converted to Sympy expression for reaction", color="red", style="normal")

    elif isinstance(e, RuntimeWarning):
        error_printer(f"\nCaught an error: ", e)
    
    elif print_trace:
        tb = traceback.extract_tb(sys.exc_info()[2])
        file_path, lineno, _, _ = tb[-1]
        file_name = os.path.basename(file_path)
        
        error_printer("Error occurred in file: ", file_name, u_end="", error_color='yellow')
        error_printer(", line: ", lineno, error_color="yellow")
        error_printer("Unexpected Error: ", e)
        error_printer("Error type: ", sys.exc_info()[0].__name__)
        message_printer("Unable to complete the query!", color="red", style="normal")