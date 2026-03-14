import logging
import sys


def setup_logging(level):
    """Handle format and output for python logger"""
    # Add root logger
    root_logger = logging.getLogger()
    root_logger.handlers.clear()
    root_logger.setLevel(level)

    # Adapt levels to C++ spdlog
    logging.addLevelName(logging.DEBUG, "debug")
    logging.addLevelName(logging.INFO, "info")
    logging.addLevelName(logging.WARNING, "warning")
    logging.addLevelName(logging.ERROR, "error")
    logging.addLevelName(logging.CRITICAL, "critical")
    
    # Terminal output
    handler = logging.StreamHandler(sys.stdout)
    terminal_formatter = logging.Formatter(
        "[%(asctime)s.%(msecs)03d] [%(name)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler.setFormatter(terminal_formatter)
    handler.setLevel(level)
    root_logger.addHandler(handler)
    
    # File output if level = DEBUG
    if level == logging.DEBUG:
        # Add manual flushhandler to ensure correct logging order with C++
        class FlushHandler(logging.FileHandler):
            def emit(self, record):
                super().emit(record)
                self.flush()

        # Set format
        fh = FlushHandler("debug.log", mode="a", encoding="utf-8", delay=True)
        file_formatter = logging.Formatter(
            "[%(asctime)s.%(msecs)03d] [%(name)s] [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        fh.setFormatter(file_formatter)
        fh.setLevel(level)
        root_logger.addHandler(fh)
    