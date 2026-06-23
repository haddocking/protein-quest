import pytest
from yarl import URL

from protein_quest.converter import Percentage, PositiveInt, Ratio, converter

items = [
    ("https://example.com", URL, URL("https://example.com")),
    (42.0, Percentage, 42.0),
    (0.42, Ratio, 0.42),
    (42, PositiveInt, 42),
]


@pytest.mark.parametrize("raw,typ,expected", items)
def test_happy_structure(raw: object, typ: type, expected: object):
    res = converter.structure(raw, typ)
    assert res == expected


@pytest.mark.parametrize("expected,typ,raw", items)
def test_happy_unstructure(expected: object, typ: type, raw: object):
    res = converter.unstructure(raw)
    assert res == expected
