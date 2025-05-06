from colorama import Fore, Back, Style, init

init( autoreset=True )

color_map = {
    "white": Fore.WHITE,
    "blue": Fore.BLUE,
    "red": Fore.RED,
    "green": Fore.GREEN,
    "yellow": Fore.YELLOW,
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