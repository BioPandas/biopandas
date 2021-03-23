import logging
import shutil
import tempfile
import sys

from os.path import isfile, join
from mkdocs.commands.build import build
from mkdocs.config import load_config

log = logging.getLogger(__name__)


def _init_asyncio_patch():
    """
    Select compatible event loop for Tornado 5+.

    As of Python 3.8, the default event loop on Windows is `proactor`,
    however Tornado requires the old default "selector" event loop.
    As Tornado has decided to leave this to users to set, MkDocs needs
    to set it. See https://github.com/tornadoweb/tornado/issues/2608.
    """
    if sys.platform.startswith("win") and sys.version_info >= (3, 8):
        import asyncio
        try:
            from asyncio import WindowsSelectorEventLoopPolicy
        except ImportError:
            pass  # Can't assign a policy which doesn't exist.
        else:
            if not isinstance(asyncio.get_event_loop_policy(), WindowsSelectorEventLoopPolicy):
                asyncio.set_event_loop_policy(WindowsSelectorEventLoopPolicy())


def _get_handler(site_dir, StaticFileHandler):

    from tornado.template import Loader

    class WebHandler(StaticFileHandler):

        def write_error(self, status_code, **kwargs):

            if status_code in (404, 500):
                error_page = '{}.html'.format(status_code)
                if isfile(join(site_dir, error_page)):
                    self.write(Loader(site_dir).load(error_page).generate())
                else:
                    super().write_error(status_code, **kwargs)

    return WebHandler


def _livereload(host, port, config, builder, site_dir, watch_theme, delay):

    # We are importing here for anyone that has issues with livereload. Even if
    # this fails, the --no-livereload alternative should still work.
    _init_asyncio_patch()
    from livereload import Server
    import livereload.handlers

    class LiveReloadServer(Server):

        def get_web_handlers(self, script):
            handlers = super().get_web_handlers(script)
            # replace livereload handler
            return [(handlers[0][0], _get_handler(site_dir, livereload.handlers.StaticFileHandler), handlers[0][2],)]

    server = LiveReloadServer()

    # Watch the documentation files, the config file and the theme files.
    server.watch(config['docs_dir'], builder, delay=delay)
    server.watch(config['config_file_path'], builder, delay=delay)

    if watch_theme:
        for d in config['theme'].dirs:
            server.watch(d, builder, delay=delay)

    # Run `serve` plugin events.
    server = config['plugins'].run_event('serve', server, config=config, builder=builder)

    server.serve(root=site_dir, host=host, port=port, restart_delay=delay)


def _static_server(host, port, site_dir):

    # Importing here to separate the code paths from the --livereload
    # alternative.
    _init_asyncio_patch()
    from tornado import ioloop
    from tornado import web

    application = web.Application([
        (r"/(.*)", _get_handler(site_dir, web.StaticFileHandler), {
            "path": site_dir,
            "default_filename": "index.html"
        }),
    ])
    application.listen(port=port, address=host)

    log.info('Running at: http://%s:%s/', host, port)
    log.info('Hold ctrl+c to quit.')
    try:
        ioloop.IOLoop.instance().start()
    except KeyboardInterrupt:
        log.info('Stopping server...')


def serve(config_file=None, dev_addr=None, strict=None, theme=None,
          theme_dir=None, livereload='livereload', watch_theme=False, wait=0, **kwargs):
    """
    Start the MkDocs development server

    By default it will serve the documentation on http://localhost:8000/ and
    it will rebuild the documentation and refresh the page automatically
    whenever a file is edited.
    """

    # Create a temporary build directory, and set some options to serve it
    # PY2 returns a byte string by default. The Unicode prefix ensures a Unicode
    # string is returned. And it makes MkDocs temp dirs easier to identify.
    site_dir = tempfile.mkdtemp(prefix='mkdocs_')

    def builder():
        log.info("Building documentation...")
        config = load_config(
            config_file=config_file,
            dev_addr=dev_addr,
            strict=strict,
            theme=theme,
            theme_dir=theme_dir,
            site_dir=site_dir,
            **kwargs
        )
        # Override a few config settings after validation
        config['site_url'] = 'http://{}/'.format(config['dev_addr'])

        live_server = livereload in ['dirty', 'livereload']
        dirty = livereload == 'dirty'
        build(config, live_server=live_server, dirty=dirty)
        return config

    try:
        # Perform the initial build
        config = builder()

        host, port = config['dev_addr']

        if livereload in ['livereload', 'dirty']:
            _livereload(host, port, config, builder, site_dir, watch_theme, delay=wait)
        else:
            _static_server(host, port, site_dir)
    finally:
        shutil.rmtree(site_dir)
