[English](README.md) | **Deutsch**

# 📈 Funktionsplotter 🚀

Ein vielseitiges Plot-Werkzeug, gebaut mit **Python**, **Tkinter** und **Matplotlib**. Von einfachen 2D-Kurven und
Parameterkurven bis zu Bode-Diagrammen und 3D-Flächen macht es aus mathematischen Ausdrücken übersichtliche,
interaktive Plots – mit automatischer Kurvenanalyse und Unterstützung für Messwerte.

👥 Teamprojekt von **Firas Mannai** und **[Yahia Bannour](https://github.com/yahiabannour)**.

| 📉 Funktionen, Kurvenanalyse und Messwerte | 📶 Bode-Diagramm | 🏔️ 3D-Fläche und Höhenlinien |
|---|---|---|
| ![2D-Plot](docs/plot-2d.png) | ![Bode-Diagramm](docs/plot-bode.png) | ![3D-Plot](docs/plot-3d.png) |

---

## ✨ Funktionen

### 🖼️ Vier Plot-Modi
* **2D-Plots**: Mehrere Funktionen gleichzeitig zeichnen, jede mit eigener Farbe, Linienbreite und Linienstil.
* **🌀 Parameterkurven**: Kurven aus x(t) und y(t) darstellen.
* **📶 Bode-Diagramme**: Betrag und Phase einer Übertragungsfunktion aus Zähler- und Nennerkoeffizienten anzeigen.
* **🏔️ 3D-Flächen & Höhenlinien**: Funktionen f(x, y) als Fläche und als Höhenlinienkarte nebeneinander zeigen.

---

### 🔬 Kurvenanalyse
* **Automatische Markierungen**: **Nullstellen**, **Extrema** und **Wendepunkte** finden und markieren –
  Nullstellen werden durch lineare Interpolation verfeinert.
* **Ableitungen**: Die numerische **1. und 2. Ableitung** jeder Funktion einblenden.

---

### 📊 Daten & Konfiguration
* **Messwerte**: x/y-Werte aus einer `.txt`- oder `.csv`-Datei laden und als Linie oder Punktwolke mit den
  Funktionen vergleichen.
* **💾 Speichern & Laden**: Alle Einstellungen in einer lesbaren Konfigurationsdatei sichern und später
  wiederherstellen.

---

### 🖱️ Interaktiv & sicher
* **Interaktive Oberfläche**: Cursor-Koordinaten live, Punkte per Klick beschriften, Zoomen und Verschieben,
  Export als Bild.
* **Modus ohne Oberfläche**: Einen Plot direkt aus einer Konfigurationsdatei erzeugen, ohne die Oberfläche zu
  öffnen.
* **🔒 Sichere Auswertung**: Ausdrücke dürfen nur freigegebene mathematische Funktionen nutzen – kein Zugriff auf
  Python-Built-ins.

---

## 🚀 Erste Schritte

### 📋 Voraussetzungen

Python 3 mit NumPy und Matplotlib:

```bash
pip install numpy matplotlib
```

### 🏃 App starten

1. **Oberfläche** (Standard):
   ```bash
   python app.py
   ```
   Öffnet das interaktive Fenster mit allen Einstellungen.

2. **Ohne Oberfläche**:
   ```bash
   python app.py config_bode.txt
   ```
   Öffnet direkt den Plot für die angegebene Konfigurationsdatei.

---

## 🛠️ Konfigurationsformat

Konfigurationsdateien nutzen ein einfaches `Schlüssel: Wert`-Format:

```text
# Zwei Funktionen mit Kurvenanalyse und Messwerten
functions: [exp(-0.5*x)*sin(x), cos(x)]
xstart: 0
xend: 10
xstep: 0.1
find_roots: True
find_extrema: True
data_file: messwerte.txt
title: Function Plot Example
```

Weitere Beispiele: [`config_bode.txt`](config_bode.txt), [`config_3d.txt`](config_3d.txt) und [`saved_config.txt`](saved_config.txt).

---

## 📂 Projektstruktur

* `app.py`: Einstiegspunkt – startet die Oberfläche, oder den Modus ohne Oberfläche, wenn eine
  Konfigurationsdatei übergeben wird.
* `gui.py`: Die Tkinter-Oberfläche mit Bedienelementen für jeden Modus und eingebetteter Matplotlib-Zeichenfläche.
* `plotter.py`: Alle Plots und die numerische Analyse, unabhängig von der Oberfläche.
* `config_parser.py`: Liest Konfigurationsdateien und wandelt die Werte in Wahrheitswerte, Zahlen und Listen um.
* `widgets.py`: Ein scrollbarer Rahmen für das lange Bedienfeld.
* `config_bode.txt`, `config_3d.txt`, `saved_config.txt`, `messwerte.txt`: Beispielkonfigurationen und Messwerte.

---

## 🧩 Unterstützte Funktionen

In den Eingabefeldern gilt die übliche mathematische Schreibweise:

* `sin`, `cos`, `tan`, `exp`, `log`, `log10`, `sqrt`, `abs`
* Konstanten: `pi`, `e`
* Potenzoperator: `^` (z. B. `x^2`)
