"""Logger that can be used together with multiprocessing funcions."""
import logging as lg
import threading
import multiprocessing as mp


class MultiprocessHandler(lg.Handler):
    """Logger class for multiprocessing functions.

    Inspired by Matt Gathu's blog post
    https://mattgathu.github.io/multiprocessing-logging-in-python/
    """

    def __init__(self, fname):
        """Multiprocessing logger."""
        lg.Handler.__init__(self)

        self._handler = lg.FileHandler(fname)
        self.queue = mp.Queue(-1)

        thread = threading.Thread(target=self.receive)
        thread.daemon = True
        thread.start()

    def setFormatter(self, fmt):
        """Set logging format."""
        lg.Handler.setFormatter(self, fmt)
        self._handler.setFormatter(fmt)

    def receive(self):
        """Receive message from queue."""
        while True:
            try:
                record = self.queue.get()
                self._handler.emit(record)
            except (KeyboardInterrupt, SystemExit):
                raise
            except EOFError:
                break
            except Exception as exc:
                print("An unexpected exception occured during logging.")
                print(str(exc))

    def send(self, msg):
        """Send message."""
        self.queue.put_nowait(msg)

    def emit(self, record):
        """Emit record."""
        try:
            self.send(record)
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception:
            self.handleError(record)

    def close(self):
        """Close file handler."""
        self._handler.close()
        lg.Handler.close(self)
