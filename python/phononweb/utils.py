import numpy as np

def open_file_phononwebsite(filename,port=8000,
                            website="http://henriquemiranda.github.io/phononwebsite",
                            host="localhost"):
    """
    take a file, detect the type and open it on the phonon website
    """
    import webbrowser
    import signal
    import sys

    #python 2 and 3 http server
    try:
        from http.server import HTTPServer, SimpleHTTPRequestHandler
    except ImportError:
        from BaseHTTPServer import HTTPServer
        # python 2 requires internal implementation
        from SimpleHTTPServer import SimpleHTTPRequestHandler

    if ".json" in filename:
        filetype = "json"
    elif ".yaml" in filename:
        filetype = "yaml"
    else: 
        filetype = "rest"

    # Add CORS header to the website
    class CORSRequestHandler (SimpleHTTPRequestHandler):
        def end_headers (self):
            self.send_header('Access-Control-Allow-Origin', "https://henriquemiranda.github.io")
            SimpleHTTPRequestHandler.end_headers(self)
        def log_message(self, format, *args):
            return

    # Quit application when SIGINT is received
    def signal_handler(signal, frame):
        sys.exit(0)

    #initialize http server thread
    print('Starting HTTP server at port %d...' % port)
    try:
        server = HTTPServer(('', port), CORSRequestHandler)
    except OSError:
        server = HTTPServer(('', port+1), CORSRequestHandler)

    server.url = 'http://{}:{}'.format(host,server.server_port)
    from threading import Thread
    t = Thread(target=server.serve_forever, daemon=True)
    t.start()

    #open website with the file
    url_filename = 'http://{}:{}/{}'.format(host,server.server_port,filename)
    url = '%s/phonon.html?%s=%s'%(website,filetype,url_filename)
    webbrowser.get(None).open_new_tab(url)

    print('Press Ctrl+C to terminate HTTP server')
    signal.signal(signal.SIGINT, signal_handler)
    signal.pause()

def estimate_band_connection(prev_eigvecs, eigvecs, prev_band_order):
    """ 
    A function to order the phonon eigenvectors taken from phonopy
    Uses the Hungarian algorithm to find the strict optimal 1-to-1 permutation
    that maximizes the overall overlap between previous and current eigenvectors.
    """
    from scipy.optimize import linear_sum_assignment
    
    # Calculate overlap metric matrix: metric[i, j] = |<prev_i | current_j>|
    metric = np.abs(np.dot(prev_eigvecs.conjugate().T, eigvecs))
    
    # linear_sum_assignment minimizes the cost. We want to maximize metric,
    # so we supply -metric to the solver.
    row_ind, col_ind = linear_sum_assignment(-metric)
    
    # row_ind is guaranteed to be [0, 1, ..., n-1]
    # col_ind[i] gives the index in the current eigvecs that corresponds to prev_eigvecs[i]
    connection_order = col_ind.tolist()

    band_order = [connection_order[x] for x in prev_band_order]
    return band_order