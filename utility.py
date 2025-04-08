from colorama import Fore, Back, Style, init

init( autoreset=True )

color_map = {
    "white": Fore.WHITE,
    "blue": Fore.BLUE,
    "red": Fore.RED,
    "green": Fore.GREEN,
    "yellow": Fore.YELLOW
}

style_map = {
    "normal": Style.NORMAL,
    "bold": Style.BRIGHT,
    "dim": Style.DIM
}


def printer(descrip, text_to_print, descript_color = "white", text_color = "blue", text_style = "normal"):

    description = f"{color_map.get(descript_color, Fore.WHITE)}{descrip}"
    text = f"{style_map.get(text_style, Style.NORMAL)}{color_map.get(text_color, Fore.BLUE)}{text_to_print}"

    print(f"{description} {text}")