# haploblock_pipeline/utils/logging.py
import logging

def setup_logger(name=None, level=logging.INFO):
    """
    Set up a logger with standard formatting.
    
    Args:
        name (str, optional): Logger name. Defaults to None (root logger).
        level (int, optional): Logging level. Defaults to logging.INFO.
    
    Returns:
        logging.Logger: Configured logger instance.
    """
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=level,
    )
    return logging.getLogger(name)

