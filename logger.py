import logging

def setup_logger():
    logging.basicConfig(
        filename="tokenizer.log",
        filemode = "a",
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s"
    )
