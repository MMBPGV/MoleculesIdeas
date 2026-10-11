import re

import pytest

from chemlabor import texte
from chemlabor.texte import TEXTE, t, t_oder, zahl


def _platzhalter(text):
    return set(re.findall(r"\{(\w+)", text))


def test_beide_sprachen_haben_dieselben_schluessel():

    nur_de = set(TEXTE["de"]) - set(TEXTE["en"])
    nur_en = set(TEXTE["en"]) - set(TEXTE["de"])

    assert not nur_de, f"Missing in en: {sorted(nur_de)}"
    assert not nur_en, f"Missing in de: {sorted(nur_en)}"


def test_platzhalter_sind_in_beiden_sprachen_gleich():

    for schluessel, de_text in TEXTE["de"].items():

        assert _platzhalter(de_text) == _platzhalter(TEXTE["en"][schluessel]), \
            f"Placeholders differ for {schluessel}"


def test_unbekannter_schluessel_liefert_schluessel():

    assert t("gibt.es.nicht") == "gibt.es.nicht"
    assert t_oder("gibt.es.nicht", "Standard") == "Standard"


def test_sprache_umschalten_und_zurueck():

    vorher = texte.sprache()

    try:

        texte.setze_sprache("de")
        assert t("phase.Automatik") == "Automatik"
        assert zahl(12345) == "12.345"

        texte.setze_sprache("en")
        assert t("phase.Automatik") == "Automatic"
        assert zahl(12345) == "12,345"

    finally:

        texte.setze_sprache(vorher)


def test_unbekannte_sprache_wird_abgelehnt():

    with pytest.raises(ValueError):
        texte.setze_sprache("xx")