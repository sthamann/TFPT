# Prime-Front-Erklärvideo (deutsch, 2026-08)

HyperFrames-Projekt für das deutschsprachige Erklärvideo der Prime Front
(2:39 min, 1920×1080, H.264 + AAC). Eingebunden auf der Website unter
`/prime-front` (Dateien in `website/public/prime-front/`).

## Inhalt (9 Frames, ~159 s)

1. Die Musik der Primzahlen (RH, kritische Linie) — 0:00
2. Der TFPT-Zugang: Geometrie zuerst (E₈ → Θ=E₄ → Hecke → ζ, v625) — 0:17
3. Die Fenstermatrix ist ein Weil-Operator (W1-Theorem, v643, Erratum) — 0:31
4. Detektor & Falsifikator (Ramanujan besteht, Epstein bricht 0,803; Matched Filter) — 0:46
5. Zwei Flächen-Sätze (det S > 0 auf 67/67, ×634; T-B 60/70) — 1:10
6. Ihara-Blaupause: der fehlende Motor Z1 (Hilbert–Pólya) — 1:29
7. Die Geometrie liefert das Maß (Λ_geo = Λ; Γ-Fluss erzwingt Massen 0,11 %) — 1:48
8. Der Korridor und die 0,53 (Levinson 0,14 %; „Erkläre die Selektion") — 2:04
9. Ehrlich: kein Beweis — ≈700 Module, 43 Lean-Theoreme — 2:24

Quelle: `big_picture_2026-08-02_de.tex` (Bogen D–G, Erratum/Theorem,
Kletter-Runde, Türklinke, die zwei 3.-August-Kapitel). Kein RH-Claim —
die Ehrlichkeits-Disziplin des Programms gilt auch im Video.

## Pipeline (faceless-explainer-Skill)

- `STORYBOARD.md` — 9 Frames mit zeitcodierten Shot-Sequenzen (VO-getaktet)
- `SCRIPT.md` — deutsche Narration (Quelle für TTS)
- `frame.md` — Design-System: broadside-Preset, remixt auf die Website-Farben
  (Slate #020617, Text #e2e8f0, Akzent Sky #38bdf8; Newsreader/Inter/JetBrains Mono)
- `scripts/german_tts.py` — deutsche TTS via **edge-tts** (de-DE-ConradNeural,
  Microsoft Neural, mit Wort-Timestamps). Grund: HeyGen offline, Kokoro hat
  kein Deutsch. Venv: `.venv-tts/`
- Captions (deutsch, burned-in) aus den Wort-Timestamps; Skin in
  `.hyperframes/caption-skin.html` (lowercase-Zwang entfernt)
- `renders/video.mp4` (27,5 MB), `renders/poster.jpg`,
  `renders/prime-front-erklaert.de.vtt` (34 Cues, satzweise)

## Neu bauen

```bash
.venv-tts/bin/python scripts/german_tts.py     # TTS + audio_meta.json
node <skill>/scripts/audio.mjs sync-durations --audio-meta ./audio_meta.json --storyboard ./STORYBOARD.md
node <skill>/scripts/captions.mjs build --storyboard ./STORYBOARD.md --audio-meta ./audio_meta.json --hyperframes . --out ./caption_groups.json
node <skill>/scripts/assemble-index.mjs --storyboard ./STORYBOARD.md --hyperframes .
node <skill>/scripts/transitions.mjs inject --storyboard ./STORYBOARD.md --hyperframes .
npx hyperframes lint && npx hyperframes validate && npx hyperframes inspect
npx hyperframes render --quality high --output renders/video.mp4
```

`<skill>` = `~/.claude/skills/faceless-explainer`. Vorschau: `npx hyperframes preview`.
