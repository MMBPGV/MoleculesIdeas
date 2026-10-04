# ==========================================================
# theme.py
# ==========================================================
#
# Farben und ttk-Stil der Oberflaeche. Farben und konfiguriere_style()
# stammen unveraendert aus der bisherigen GUI (app.py); ergaenze_style()
# kommt neu dazu (Combobox, Menubutton, Scrollbar).

from tkinter import ttk


FARBEN = {
    "bg": "#1e1e1e",
    "panel": "#252526",
    "panel_hell": "#2d2d30",
    "rand": "#3f3f46",
    "text": "#e0e0e0",
    "text_gedimmt": "#9a9a9a",
    "akzent": "#4fc3f7",
    "akzent_dunkel": "#2f9ee0",
    "gefahr": "#7a3030",
    "gefahr_hell": "#a23f3f",
    "auswahl": "#3a5f8a",
}

ELEMENT_FARBEN = {
    "H": "#e8e8e8", "C": "#4a4a4a", "N": "#3060c8", "O": "#c83030",
    "F": "#3fae3f", "Cl": "#2f9e5f", "P": "#e08820", "S": "#d8c020",
    "Si": "#a8763f",
    # in der bisherigen Tabelle fehlten diese drei - sie erschienen in der Standardfarbe
    "B": "#e8a0a0",
    "Br": "#a03030",
    "I": "#8030a0",
}


def konfiguriere_style(root):

    f = FARBEN

    root.configure(bg=f["bg"])

    style = ttk.Style(root)
    style.theme_use("clam")

    style.configure(
        ".",
        background=f["bg"],
        foreground=f["text"],
        fieldbackground=f["panel_hell"],
        bordercolor=f["rand"],
        font=("Segoe UI", 10)
    )

    style.configure("TFrame", background=f["bg"])

    style.configure("TLabel", background=f["bg"], foreground=f["text"])

    style.configure(
        "TLabelframe",
        background=f["bg"],
        foreground=f["text"],
        bordercolor=f["rand"]
    )

    style.configure(
        "TLabelframe.Label",
        background=f["bg"],
        foreground=f["akzent"],
        font=("Segoe UI", 10, "bold")
    )

    style.configure(
        "TButton",
        background=f["panel_hell"],
        foreground=f["text"],
        bordercolor=f["rand"],
        padding=6
    )

    style.map(
        "TButton",
        background=[("active", f["auswahl"]), ("disabled", f["panel"])],
        foreground=[("disabled", f["text_gedimmt"])]
    )

    # Hervorgehobene Buttons fuer die wichtigsten Aktionen (Start)
    # bzw. eine potenziell folgenreiche Aktion (Stoppen) - macht auf
    # einen Blick klar, welcher Button "los geht's" und welcher
    # "haelt an" bedeutet, statt dass beide gleich aussehen.
    style.configure(
        "Accent.TButton",
        background=f["akzent_dunkel"],
        foreground="#ffffff",
        bordercolor=f["akzent_dunkel"],
        padding=6,
        font=("Segoe UI", 10, "bold")
    )

    style.map(
        "Accent.TButton",
        background=[("active", f["akzent"]), ("disabled", f["panel"])]
    )

    style.configure(
        "Danger.TButton",
        background=f["gefahr"],
        foreground="#ffffff",
        bordercolor=f["gefahr"],
        padding=6
    )

    style.map(
        "Danger.TButton",
        background=[("active", f["gefahr_hell"]), ("disabled", f["panel"])]
    )

    style.configure(
        "TEntry",
        fieldbackground=f["panel_hell"],
        foreground=f["text"],
        insertcolor=f["text"],
        bordercolor=f["rand"]
    )

    style.configure("TCheckbutton", background=f["bg"], foreground=f["text"])

    style.map(
        "TCheckbutton",
        background=[("active", f["bg"])],
        foreground=[("disabled", f["text_gedimmt"])]
    )

    style.configure(
        "TNotebook",
        background=f["bg"],
        bordercolor=f["rand"]
    )

    style.configure(
        "TNotebook.Tab",
        background=f["panel"],
        foreground=f["text_gedimmt"],
        padding=(14, 7)
    )

    style.map(
        "TNotebook.Tab",
        background=[("selected", f["panel_hell"])],
        foreground=[("selected", f["akzent"])]
    )

    style.configure(
        "Treeview",
        background=f["panel_hell"],
        fieldbackground=f["panel_hell"],
        foreground=f["text"],
        bordercolor=f["rand"],
        borderwidth=0,
        rowheight=24
    )

    style.configure(
        "Treeview.Heading",
        background=f["panel"],
        foreground=f["akzent"],
        bordercolor=f["rand"],
        font=("Segoe UI", 9, "bold")
    )

    style.map(
        "Treeview",
        background=[("selected", f["auswahl"])],
        foreground=[("selected", "#ffffff")]
    )

    for orientierung in ("Vertical", "Horizontal"):

        style.configure(
            f"{orientierung}.TScrollbar",
            background=f["panel"],
            troughcolor=f["bg"],
            bordercolor=f["rand"],
            arrowcolor=f["text"]
        )


def ergaenze_style(root):
    """Dunkler Stil fuer Widgets, die die bisherige GUI nicht benutzte."""

    style = ttk.Style(root)

    style.configure(
        "TCombobox",
        fieldbackground=FARBEN["panel_hell"],
        background=FARBEN["panel_hell"],
        foreground=FARBEN["text"],
        arrowcolor=FARBEN["text"],
        bordercolor=FARBEN["rand"],
    )

    style.map(
        "TCombobox",
        fieldbackground=[("readonly", FARBEN["panel_hell"])],
        foreground=[("readonly", FARBEN["text"])],
        selectbackground=[("readonly", FARBEN["panel_hell"])],
        selectforeground=[("readonly", FARBEN["text"])],
    )

    root.option_add("*TCombobox*Listbox.background", FARBEN["panel_hell"])
    root.option_add("*TCombobox*Listbox.foreground", FARBEN["text"])
    root.option_add("*TCombobox*Listbox.selectBackground", FARBEN["auswahl"])

    style.configure(
        "TMenubutton",
        background=FARBEN["panel_hell"],
        foreground=FARBEN["text"],
        bordercolor=FARBEN["rand"],
        padding=(10, 5),
    )

    style.map("TMenubutton", background=[("active", FARBEN["auswahl"])])

    style.configure(
        "Vertical.TScrollbar",
        background=FARBEN["panel_hell"],
        troughcolor=FARBEN["panel"],
        bordercolor=FARBEN["panel"],
        arrowcolor=FARBEN["text"],
    )

    style.configure("Gedimmt.TLabel", foreground=FARBEN["text_gedimmt"])
    style.configure("Kennzahl.TLabel", font=("", 16, "bold"), foreground=FARBEN["akzent"])
