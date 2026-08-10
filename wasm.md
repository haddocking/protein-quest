# wasm version of protein-quest

TODO remove once PR is ready to merge

Via mkdocs plugin
```shell
uv add --group docs mkdocs-jupyterlite
uv run mkdocs build
caddy file-server --root site --listen 0.0.0.0:8080
```

Or build own lab env and use iframe in mkdocs page

```shell
jupyter lite build --debug --contents files --piplite-wheels /tmp/tmpv9apl6j9/wheels/wheels/gemmi-0.7.6.dev0-cp314-cp314-pyemscripten_2026_0_wasm32.whl
           --piplite-wheels /tmp/tmpv9apl6j9/wheels/wheels/mmcif-1.1.1-cp314-cp314-pyemscripten_2026_0_wasm32.whl --piplite-wheels /tmp/tmpv9apl6j9/wheels/wheels/protein_quest-1.6.0-py3-none-any.whl
           --no-libarchive --apps notebooks --no-unused-shared-packages --output-dir /tmp/tmp1h02_xi7


uv build --wheel --out-dir docs/notebooks/
rm -rf scratch/nbwasm && uvx --from  jupyterlite-core --with jupyter-server --with jupyterlite-pyodide-kernel jupyter lite build --with libarchive-c --contents docs/notebooks/ --output-dir scratch/nbwasm -apps lab --no-unused-shared-packages

# rm -rf scratch/nbwasm && uvx --from  jupyterlite-core --with jupyter-server --with jupyterlite-pyodide-kernel jupyter lite build \
# --with libarchive-c --contents docs/notebooks/ --output-dir scratch/nbwasm -apps lab --no-unused-shared-packages \
# --piplite-wheels docs/notebooks/gemmi-0.7.6.dev0-cp314-cp314-pyemscripten_2026_0_wasm32.whl \
# --piplite-wheels docs/notebooks/mmcif-1.1.1-cp314-cp314-pyemscripten_2026_0_wasm32.whl \
# --piplite-wheels docs/notebooks/protein_quest-1.6.0-py3-none-any.whl
# does not make `%pip install protein-quest` work at prefers pypi

caddy file-server --root scratch/nbwasm --listen 0.0.0.0:8081
```

```python
import micropip

gemmi_whl = "http://localhost:8080/gemmi-0.7.6.dev0-cp314-cp314-pyemscripten_2026_0_wasm32.whl"
await micropip.install(gemmi_whl)
mmcif_whl = "http://localhost:8080/mmcif-1.1.1-cp314-cp314-pyemscripten_2026_0_wasm32.whl"
await micropip.install(mmcif_whl)
pq_whl = "http://localhost:8080/protein_quest-1.6.0-py3-none-any.whl"
await micropip.install(pq_whl)
```


## psutil

Make non-dep in wasm and added workaround

```shell
uv build --wheel
cp dist/protein_quest-1.6.0-py3-none-any.whl docs/notebookswasm/
```

## gemmi

Make wasm using cibuildwheel for python 3.14.2

```shell
gh extension install https://github.com/nektos/gh-act
cp .github/workflows/wheels2.yml .github/workflows/wheels_wasm.yml
# pick medium sized runner image
gh act workflow_dispatch -j build_wasm_wheels --artifact-server-path $PWD/.artifacts
unzip -d ../protein-quest/docs/notebookswasm/ .artifacts/1/cibw-wheels--/cibw-wheels--.zip
```

`.github/workflows/wheels2wasm.yml`:

```yaml
# Use cibuildwheel v2 to build wheels for Python 3.
# Based on https://cibuildwheel.readthedocs.io/en/stable/setup/

name: Wheels2 wasm only

on:
  workflow_dispatch:

jobs:
  build_wasm_wheels:
    name: Wheels wasm
    runs-on: ubuntu-latest

    steps:
    - uses: actions/checkout@v5

    - name: Build wheels
      uses: pypa/cibuildwheel@v4.2.0
      env:
        CIBW_BUILD: "cp314-pyodide_wasm32"
        CIBW_PLATFORM: pyodide
        CIBW_ARCHS: wasm32
        # Skip tests, as some others fail for pyodide
        CIBW_TEST_COMMAND_PYODIDE: "true"

    - uses: actions/upload-artifact@v4
      with:
        name: cibw-wheels-wasm
        path: ./wheelhouse/*.whl
```

test
```python
gemmi_whl = "http://localhost:8080/gemmi-0.7.6.dev0-cp314-cp314-pyemscripten_2026_0_wasm32.whl"
await micropip.install(gemmi_whl)
import urllib.request
import gemmi

url = "https://www.ebi.ac.uk/pdbe/entry-files/download/3jrs_updated.cif"
structure = gemmi.cif.read_string(urllib.request.urlopen(url).read())
structure
```

## mmcif

```shell
git clone  --recurse-submodules https://github.com/rcsb/py-mmcif
cd py-mmcif
mkdir -p .github/workflows
cp ../gemmi/.github/workflows/wheels2wasm.yml .github/workflows/wheels_wasm.yml
rm  ~/.config/act/actrc 
# use large image with cmake
gh act workflow_dispatch -j build_wasm_wheels --artifact-server-path $PWD/.artifacts
unzip -d ../protein-quest/docs/notebookswasm/ .artifacts/1/cibw-wheels-wasm/cibw-wheels-wasm.zip
```

Needed changes to code see https://github.com/i-VRESSE/py-mmcif fork.

test
```
from mmcif.api.DictionaryApi import DictionaryApi
from mmcif.io.BinaryCifReader import BinaryCifReader
from mmcif.io.BinaryCifWriter import BinaryCifWriter
from mmcif.io.PdbxReader import PdbxReader
from mmcif.io.PdbxWriter import PdbxWriter
```

## distributed

Made map_with_progress work when distributed is available or not.


## Uniprot notebook

```python
---------------------------------------------------------------------------
RuntimeError                              Traceback (most recent call last)
Cell In[5], line 1
----> 1 uniprot_accessions = search4uniprot(query, limit=200)

File /lib/python3.14/site-packages/protein_quest/uniprot.py:342, in search4uniprot(query, limit, timeout)
    339 sparql_query = _build_sparql_query_uniprot(query, limit)
    340 logger.info("Executing SPARQL query for UniProt: %s", sparql_query)
--> 342 raw_results = execute_sparql_search(
    343     sparql_query=sparql_query,
    344     timeout=timeout,
    345 )
    346 _limit_check("Search for uniprot accessions", limit, len(raw_results))
    347 return {result["protein"]["value"].split("/")[-1] for result in raw_results}

File /lib/python3.14/site-packages/protein_quest/sparql.py:110, in execute_sparql_search(sparql_query, timeout)
    107     sparql.setMethod("POST")
    109 sparql.setQuery(sparql_query)
--> 110 rawresults = sparql.queryAndConvert()
    111 if not isinstance(rawresults, dict):
    112     msg = f"Expected rawresults to be a dict, but got {type(rawresults)}"

File /lib/python3.14/site-packages/SPARQLWrapper/Wrapper.py:967, in SPARQLWrapper.queryAndConvert(self)
    962 def queryAndConvert(self) -> "QueryResult.ConvertResult":
    963     """Macro like method: issue a query and return the converted results.
    964 
    965     :return: the converted query result. See the conversion methods for more details.
    966     """
--> 967     res = self.query()
    968     return res.convert()

File /lib/python3.14/site-packages/SPARQLWrapper/Wrapper.py:960, in SPARQLWrapper.query(self)
    942 def query(self) -> "QueryResult":
    943     """
    944     Execute the query.
    945     Exceptions can be raised if either the URI is wrong or the HTTP sends back an error (this is also the
   (...)    958     :rtype: :class:`QueryResult` instance
    959     """
--> 960     return QueryResult(self._query())

File /lib/python3.14/site-packages/SPARQLWrapper/Wrapper.py:924, in SPARQLWrapper._query(self)
    922 try:
    923     if self.timeout:
--> 924         response = urlopener(request, timeout=self.timeout)
    925     else:
    926         response = urlopener(request)

File /lib/python314.zip/urllib/request.py:187, in urlopen(url, data, timeout, context)
    183     elif _opener is None:
    184         _opener = opener = build_opener()
    185     else:
    186         opener = _opener
--> 187     return opener.open(url, data, timeout)

File /lib/python314.zip/urllib/request.py:487, in OpenerDirector.open(self, fullurl, data, timeout)
    483             meth = getattr(processor, meth_name)
    484             req = meth(req)
    485 
    486         sys.audit('urllib.Request', req.full_url, req.data, req.headers, req.get_method())
--> 487         response = self._open(req, data)
    488 
    489         # post-process response
    490         meth_name = protocol+"_response"

File /lib/python314.zip/urllib/request.py:504, in OpenerDirector._open(self, req, data)
    500         if result:
    501             return result
    502 
    503         protocol = req.type
--> 504         result = self._call_chain(self.handle_open, protocol, protocol +
    505                                   '_open', req)
    506         if result:
    507             return result

File /lib/python314.zip/urllib/request.py:464, in OpenerDirector._call_chain(self, chain, kind, meth_name, *args)
    460         # could.  Otherwise, they return the response.
    461         handlers = chain.get(kind, ())
    462         for handler in handlers:
    463             func = getattr(handler, meth_name)
--> 464             result = func(*args)
    465             if result is not None:
    466                 return result

File /lib/python314.zip/urllib/request.py:1369, in HTTPSHandler.https_open(self, req)
   1368         def https_open(self, req):
-> 1369             return self.do_open(http.client.HTTPSConnection, req,
   1370                                 context=self._context)

File /lib/python314.zip/urllib/request.py:1328, in AbstractHTTPHandler.do_open(self, http_class, req, **http_conn_args)
   1324                 raise URLError(err)
   1325             r = h.getresponse()
   1326         except:
   1327             h.close()
-> 1328             raise
   1329 
   1330         # If the server does not send us a 'Connection: close' header,
   1331         # HTTPConnection assumes the socket should be left open. Manually

File /lib/python314.zip/http/client.py:1358, in HTTPConnection.request(self, method, url, body, headers, encode_chunked)
   1355     def request(self, method, url, body=None, headers={}, *,
   1356                 encode_chunked=False):
   1357         """Send a complete request to the server."""
-> 1358         self._send_request(method, url, body, headers, encode_chunked)

File /lib/python314.zip/http/client.py:1404, in HTTPConnection._send_request(self, method, url, body, headers, encode_chunked)
   1400         if isinstance(body, str):
   1401             # RFC 2616 Section 3.7.1 says that text default has a
   1402             # default charset of iso-8859-1.
   1403             body = _encode(body, 'body')
-> 1404         self.endheaders(body, encode_chunked=encode_chunked)

File /lib/python314.zip/http/client.py:1353, in HTTPConnection.endheaders(self, message_body, encode_chunked)
   1349         if self.__state == _CS_REQ_STARTED:
   1350             self.__state = _CS_REQ_SENT
   1351         else:
   1352             raise CannotSendHeader()
-> 1353         self._send_output(message_body, encode_chunked=encode_chunked)

File /lib/python314.zip/http/client.py:1113, in HTTPConnection._send_output(self, message_body, encode_chunked)
   1109         """
   1110         self._buffer.extend((b"", b""))
   1111         msg = b"\r\n".join(self._buffer)
   1112         del self._buffer[:]
-> 1113         self.send(msg)
   1114 
   1115         if message_body is not None:
   1116 

File /lib/python314.zip/http/client.py:1057, in HTTPConnection.send(self, data)
   1053         """
   1054 
   1055         if self.sock is None:
   1056             if self.auto_open:
-> 1057                 self.connect()
   1058             else:
   1059                 raise NotConnected()
   1060 

File /lib/python314.zip/http/client.py:1499, in HTTPSConnection.connect(self)
   1495                 server_hostname = self._tunnel_host
   1496             else:
   1497                 server_hostname = self.host
   1498 
-> 1499             self.sock = self._context.wrap_socket(self.sock,
   1500                                                   server_hostname=server_hostname)

File /lib/python314.zip/ssl.py:312, in SSLContext.wrap_socket(self, sock, server_side, do_handshake_on_connect, suppress_ragged_eofs, server_hostname, session)
    308         suppress_ragged_eofs=True,
    309         server_hostname=None,
    310         session=None,
    311     ):
--> 312         return SSLSocket._create(
    313             sock,
    314             server_side=server_side,
    315             do_handshake_on_connect=do_handshake_on_connect,

File /lib/python314.zip/ssl.py:739, in SSLSocket._create(cls, sock, server_side, do_handshake_on_connect, suppress_ragged_eofs, server_hostname, context, session)
    735         ssl_sock._suppress_ragged_eofs = suppress_ragged_eofs
    736         ssl_sock._session = session
    737 
    738         if do_handshake_on_connect:
--> 739             ssl_sock.do_handshake()
    740         return ssl_sock

File /lib/python314.zip/ssl.py:842, in SSLSocket.do_handshake(self, block)
    838     def do_handshake(self, block=False):
    839         try:
    840             from pyodide_js._api import _nodeSock
    841         except ImportError:
--> 842             raise RuntimeError("TLS not supported in this environment") from None
    843 
    844         result = _nodeSock.startTls(self.fileno())
    845         if result < 0:

RuntimeError: TLS not supported in this environment
```

## PDBe notebook

```python
/lib/python3.14/site-packages/tqdm/asyncio.py:67: TqdmMonitorWarning: tqdm:disabling monitor support (monitor_interval = 0) due to:
can't start new thread
  yield from cls(asyncio.as_completed(fs, timeout=timeout, **kwargs),
Downloading PDBe mmCIF files:   0%|          | 0/4 [00:00<?, ?it/s]
---------------------------------------------------------------------------
gaierror                                  Traceback (most recent call last)
File /lib/python3.14/site-packages/aiohttp/connector.py:1574, in TCPConnector._create_direct_connection(self, req, traces, timeout, client_error)
   1570 try:
   1571     # Cancelling this lookup should not cancel the underlying lookup
   1572     #  or else the cancel event will get broadcast to all the waiters
   1573     #  across all connections.
-> 1574     hosts = await self._resolve_host(host, port, traces=traces)
   1575 except OSError as exc:

File /lib/python3.14/site-packages/aiohttp/connector.py:1190, in TCPConnector._resolve_host(self, host, port, traces)
   1189 try:
-> 1190     return await asyncio.shield(resolved_host_task)
   1191 except asyncio.CancelledError:

File /lib/python3.14/site-packages/aiohttp/connector.py:1221, in TCPConnector._resolve_host_with_throttle(self, key, host, port, futures, traces)
   1219         await trace.send_dns_resolvehost_start(host)
-> 1221 addrs = await self._resolver.resolve(host, port, family=self._family)
   1222 if traces:

File /lib/python3.14/site-packages/aiohttp/resolver.py:40, in ThreadedResolver.resolve(self, host, port, family)
     37 async def resolve(
     38     self, host: str, port: int = 0, family: socket.AddressFamily = socket.AF_INET
     39 ) -> List[ResolveResult]:
---> 40     infos = await self._loop.getaddrinfo(
     41         host,
     42         port,
     43         type=socket.SOCK_STREAM,
     44         family=family,
     45         flags=_AI_ADDRCONFIG,
     46     )
     48     hosts: List[ResolveResult] = []

File /lib/python314.zip/pyodide/webloop.py:802, in WebLoop.getaddrinfo(self, host, port, family, type, proto, flags)
    800         import socket
    801 
--> 802         return await self.run_in_executor(
    803             None, socket.getaddrinfo, host, port, family, type, proto, flags

File /lib/python314.zip/pyodide/webloop.py:530, in WebLoop.run_in_executor(self, executor, func, *args)
    528         except BaseException as e:
    529             fut.set_exception(e)
--> 530         return fut

File /lib/python314.zip/socket.py:983, in getaddrinfo(host, port, family, type, proto, flags)
    981     # and socket type values to enum constants.
    982     addrlist = []
--> 983     for res in _socket.getaddrinfo(host, port, family, type, proto, flags):
    984         af, socktype, proto, canonname, sa = res

gaierror: [Errno -2] Name does not resolve

The above exception was the direct cause of the following exception:

ClientConnectorDNSError                   Traceback (most recent call last)
Cell In[5], line 1
----> 1 files = await fetch(ids, save_dir)
      2 files

File /lib/python3.14/site-packages/protein_quest/pdbe/fetch.py:63, in fetch(ids, save_dir, archived, max_parallel_downloads, cacher)
     60 urls = list(id2urls.values())
     61 id2paths = {pdb_id: save_dir / fn for pdb_id, (_, fn) in id2urls.items()}
---> 63 await retrieve_files(urls, save_dir, max_parallel_downloads, desc="Downloading PDBe mmCIF files", cacher=cacher)
     64 return id2paths

File /lib/python3.14/site-packages/protein_quest/utils.py:310, in retrieve_files(urls, save_dir, max_parallel_downloads, retries, total_timeout, desc, cacher, chunk_size, gzip_files, raise_for_not_found)
    296 async with friendly_session(retries, total_timeout) as session:
    297     tasks = [
    298         _retrieve_file(
    299             session=session,
   (...)    308         for url, filename, *rest in urls
    309     ]
--> 310     raw_files: list[Path | None] = await tqdm.gather(*tasks, desc=desc)
    311     return [f for f in raw_files if f is not None]

File /lib/python3.14/site-packages/tqdm/asyncio.py:79, in tqdm_asyncio.gather(cls, loop, timeout, total, *fs, **tqdm_kwargs)
     76     return i, await f
     78 ifs = [wrap_awaitable(i, f) for i, f in enumerate(fs)]
---> 79 res = [await f for f in cls.as_completed(ifs, loop=loop, timeout=timeout,
     80                                          total=total, **tqdm_kwargs)]
     81 return [i for _, i in sorted(res)]

File /lib/python314.zip/asyncio/tasks.py:618, in _AsCompletedIterator._wait_for_one(self, resolve)
    614         f = await self._done.get()
    615         if f is None:
    616             # Dummy value from _handle_timeout().
    617             raise exceptions.TimeoutError
--> 618         return f.result() if resolve else f

File /lib/python3.14/site-packages/tqdm/asyncio.py:76, in tqdm_asyncio.gather.<locals>.wrap_awaitable(i, f)
     75 async def wrap_awaitable(i, f):
---> 76     return i, await f

File /lib/python3.14/site-packages/protein_quest/utils.py:362, in _retrieve_file(session, url, save_path, semaphore, cacher, chunk_size, gzip_files, raise_for_not_found)
    358 auto_decompress = not gzip_files
    359 headers = {"Accept-Encoding": "gzip"}
    360 async with (
    361     semaphore,
--> 362     session.get(url, headers=headers, auto_decompress=auto_decompress) as resp,
    363 ):
    364     if not raise_for_not_found and resp.status == 404:
    365         logger.debug(f"File not found at {url}, skipping download.")

File /lib/python3.14/site-packages/aiohttp_retry/client.py:158, in _RequestContext.__aenter__(self)
    157 async def __aenter__(self) -> ClientResponse:
--> 158     return await self._do_request()

File /lib/python3.14/site-packages/aiohttp_retry/client.py:119, in _RequestContext._do_request(self)
    116 except IndexError:
    117     params = self._params_list[-1]
--> 119 response: ClientResponse = await self._request_func(
    120     params.method,
    121     params.url,
    122     headers=params.headers,
    123     trace_request_ctx={
    124         "current_attempt": current_attempt,
    125         **(params.trace_request_ctx or {}),
    126     },
    127     **(params.kwargs or {}),
    128 )
    130 debug_message = f"Retrying after response code: {response.status}"
    131 skip_retry = await self._is_skip_retry(current_attempt, response)

File /lib/python3.14/site-packages/aiohttp/client.py:788, in ClientSession._request(self, method, str_or_url, params, data, json, cookies, headers, skip_auto_headers, auth, allow_redirects, max_redirects, compress, chunked, expect100, raise_for_status, read_until_eof, proxy, proxy_auth, timeout, verify_ssl, fingerprint, ssl_context, ssl, server_hostname, proxy_headers, trace_request_ctx, read_bufsize, auto_decompress, max_line_size, max_field_size, max_headers, middlewares)
    785     handler = _connect_and_send_request
    787 try:
--> 788     resp = await handler(req)
    789 # Client connector errors should not be retried
    790 except (
    791     ConnectionTimeoutError,
    792     ClientConnectorError,
    793     ClientConnectorCertificateError,
    794     ClientConnectorSSLError,
    795 ):

File /lib/python3.14/site-packages/aiohttp/client.py:742, in ClientSession._request.<locals>._connect_and_send_request(req)
    740 assert self._connector is not None
    741 try:
--> 742     conn = await self._connector.connect(
    743         req, traces=traces, timeout=real_timeout
    744     )
    745 except asyncio.TimeoutError as exc:
    746     raise ConnectionTimeoutError(
    747         f"Connection timeout to host {req.url}"
    748     ) from exc

File /lib/python3.14/site-packages/aiohttp/connector.py:672, in BaseConnector.connect(self, req, traces, timeout)
    670     for trace in traces:
    671         await trace.send_connection_create_start()
--> 672 proto = await self._create_connection(req, traces, timeout)
    673 if traces:
    674     for trace in traces:

File /lib/python3.14/site-packages/aiohttp/connector.py:1251, in TCPConnector._create_connection(self, req, traces, timeout)
   1249     _, proto = await self._create_proxy_connection(req, traces, timeout)
   1250 else:
-> 1251     _, proto = await self._create_direct_connection(req, traces, timeout)
   1253 return proto

File /lib/python3.14/site-packages/aiohttp/connector.py:1580, in TCPConnector._create_direct_connection(self, req, traces, timeout, client_error)
   1577         raise
   1578     # in case of proxy it is not ClientProxyConnectionError
   1579     # it is problem of resolving proxy ip itself
-> 1580     raise ClientConnectorDNSError(req.connection_key, exc) from exc
   1582 last_exc: Optional[Exception] = None
   1583 addr_infos = self._convert_hosts_to_addr_infos(hosts)

ClientConnectorDNSError: Cannot connect to host www.ebi.ac.uk:443 ssl:default [Name does not resolve]
```