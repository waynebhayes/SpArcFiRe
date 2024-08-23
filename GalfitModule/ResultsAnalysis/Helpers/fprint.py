def fprint(input_str, fill_char = "*", fill_len = 100):
    input_str = f" {input_str} "
    print()
    print(f"{input_str:{fill_char}^{fill_len}}")
    print()