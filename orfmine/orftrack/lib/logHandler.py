import logging
import logging.config


class Logger:
  def __init__(self, name, outpath=None):
    if outpath:
        config = get_config(outpath=outpath)
        logging.config.dictConfig(config)
    self.logger = logging.getLogger(name)

  def deco(self, log):
    return '-' * len(log)

  def title(self, log):
      self.logger.info('')
      self.logger.info('')
      self.logger.info('')
      deco = self.deco(log)
      # self.logger.info(deco)
      self.logger.info(log.upper())
      self.logger.info(deco)
      self.logger.info('')

  def debug(self, log):
      self.logger.debug(log)

  def info(self, log):
      self.logger.info(log)

  def warning(self, log):
      self.logger.warning(log)

  def error(self, log):
      self.logger.error(log)

  def critical(self, log):
      self.logger.critical(log)


def get_logger(name, outpath=None):
    if outpath:
        config = get_config(outpath=outpath)
        logging.config.dictConfig(config)

    return logging.getLogger(name)



def get_config(outpath='./'):
  config = {
      "version": 1,
      "disable_existing_loggers": False,
      "formatters": {
        "detailed": {
          "format": "%(asctime)-12s %(name)-24s %(levelname)s:\t\t %(message)s",
          "datefmt": "%H:%M:%S"
        },
        "simple": {
          "format": "%(message)s"
        }
      },

      "handlers": {
        "console": {
          "class": "logging.StreamHandler",
          "level": "DEBUG",
          "formatter": "simple"
        },

        "info_handler": {
          "class": "logging.FileHandler",
          "level": "INFO",
          "formatter": "simple",
          "filename": outpath+"info.log",
          "mode": "w+",
          "encoding": "utf8",
          "delay": False
        },

        "summary_handler": {
          "class": "logging.FileHandler",
          "level": "INFO",
          "formatter": "simple",
          "filename": outpath+"summary.log",
              "mode": "w+",
          "encoding": "utf8",
          "delay": False
        }
      },

      "loggers": {
        "orftrack.lib": {
          "level": "DEBUG",
          "handlers": ["console", "info_handler"],
          "propagate": False
        },
        "summary": {
          "level": "INFO",
          "handlers": ["console", "summary_handler"],
          "propagate": False
        }
      },

      "root": {
        "level": "NOTSET",
        "handlers": ["console", "info_handler"],
        "propagate": False
      }
  }

  return config
