"use client";

import { useEffect, useRef, useState } from "react";
import { Play } from "lucide-react";
import { SectionHeader } from "../SectionHeader";
import { cn, SITE_URL } from "@/lib/utils";

const VIDEO_SRC = "/prime-front/prime-front-erklaert.mp4";
const POSTER_SRC = "/prime-front/prime-front-erklaert-poster.jpg";
const CAPTIONS_SRC = "/prime-front/prime-front-erklaert.de.vtt";

/** Chapter markers — seconds match the film's nine frames (STORYBOARD.md). */
const CHAPTERS: { t: number; label: string }[] = [
  { t: 0, label: "Die Musik der Primzahlen" },
  { t: 17, label: "Geometrie zuerst (E₈ → ζ)" },
  { t: 31, label: "Das W1-Theorem" },
  { t: 46, label: "Detektor & Falsifikator" },
  { t: 70, label: "Zwei Flächen-Sätze" },
  { t: 89, label: "Die Ihara-Blaupause (Z1)" },
  { t: 108, label: "Das Maß aus der Geometrie" },
  { t: 124, label: "Der Korridor & die 0,53" },
  { t: 144, label: "Ehrlich: kein Beweis" },
];

/**
 * The German narration, grouped by chapter. Mirrored from the film's caption
 * track — rendered as visible (collapsible) text so it is accessible to screen
 * readers and indexable (the video pixels are not). Keep in sync with
 * /prime-front/prime-front-erklaert.de.vtt.
 */
const TRANSCRIPT: { heading: string; body: string }[] = [
  {
    heading: "Die Musik der Primzahlen (0:00)",
    body: "Die Primzahlen wirken wie Zufall. Doch ihre Verteilung folgt einem verborgenen Orchester: den Nullstellen der Zetafunktion. Die Riemannsche Vermutung sagt: Alle liegen auf einer einzigen Linie. Unbewiesen — seit über hundertfünfundsechzig Jahren.",
  },
  {
    heading: "Geometrie zuerst (0:17)",
    body: "Die TFPT baut aus zwei Axiomen das E-acht-Gitter. Seine Zählfunktion kennt die Primzahlen bereits: Jede wirkt als eigener Prüfkanal. Und die Zetafunktion erscheint als Schatten dieser Geometrie.",
  },
  {
    heading: "Das W1-Theorem (0:31)",
    body: "Aus dieser Buchführung entstand eine Fenstermatrix. Und sie entpuppte sich als etwas Klassisches: Wort für Wort Suzukis Weil-Operator. Das ist das W-eins-Theorem — maschinell verifiziert, mit einem ehrlich dokumentierten Erratum.",
  },
  {
    heading: "Detektor & Falsifikator (0:46)",
    body: "Dieselbe Form ist ein Detektor. Kalibriert an gelösten Welten: Ramanujan-Graphen — dort ist das Analogon bewiesen — bestehen den Test. Epstein — mit echten Nullstellen neben der Linie — bricht, exakt wie vorhergesagt. Und der Matched Filter macht es konstruktiv: Jede Off-Line-Nullstelle würde einen nachrechenbaren Zeugen erzeugen. Das Objekt ist falsifizierbar.",
  },
  {
    heading: "Zwei Flächen-Sätze (1:10)",
    body: "Zwei Sätze stehen auf der ganzen Fläche. Das Vorzeichen der Determinante: unbedingt bewiesen, auf allen siebenundsechzig Fenstern — die hundertjährige nullstellenfreie Zone versperrt den einzigen Fluchtweg. Und die Marge: sechzig von siebzig Fenstern, geschlossen mit zitierter Klassik.",
  },
  {
    heading: "Die Ihara-Blaupause (1:29)",
    body: "Dann der Befund vom dritten August: Im Graphen-Labor, wo das Analogon bewiesen ist, existiert die Ziel-Zerlegung exakt — Quadratsumme plus Defekt. Unsere Fensterform ist baugleich. Es fehlt genau ein Bauteil: der Motor. Hilbert–Pólya, in Fensterkoordinaten.",
  },
  {
    heading: "Das Maß aus der Geometrie (1:48)",
    body: "Und die Geometrie liefert schon das Maß: Die Primzahl-Atome lassen sich zirkelfrei aus reiner Gitterzählung ablesen. Der Gamma-Fluss erzwingt ihre Massen auf Promille — und ihre Positionen. Aber ehrlich: als Prüfstand, nicht als Generator.",
  },
  {
    heading: "Der Korridor und die 0,53 (2:04)",
    body: "Der Stand heute: Jede Masse lebt in einem Korridor mit exakt berechenbaren Rändern. Die Arithmetik wählt nicht den Rand — sie wählt einen inneren Punkt, bei null Komma dreiundfünfzig. Ein Energie-Extremum trifft ihn bis auf Promille. Die offene Frage: Erkläre die Selektion im Korridor.",
  },
  {
    heading: "Ehrlich: kein Beweis (2:24)",
    body: "Kein Beweis. Das sagt dieses Programm selbst, an jeder Stelle. Aber die Frage war nie so klein — und nie so präzise. Rund siebenhundert Module, jede Zahl maschinell geprüft. Das ist die Prime Front.",
  },
];

const fmt = (t: number) => `${Math.floor(t / 60)}:${String(t % 60).padStart(2, "0")}`;

const videoJsonLd = {
  "@context": "https://schema.org",
  "@type": "VideoObject",
  name: "Die Prime Front — komplett erklärt (deutsch)",
  description:
    "A ~2.5 minute German-language explainer of TFPT's prime / zeta line, end to end: the Riemann Hypothesis in one image (all zeros on one line), the window matrix born from the E8 bookkeeping and proved identical to Suzuki's Weil operator (the W1 theorem, v643), the calibrated detector / falsifier (Ramanujan passes, Epstein breaks at the predicted 0.803; matched-filter witness for any off-line zero), the two surface theorems (det S > 0 unconditionally on 67/67 windows; the T-B margin closed on 60/70 with cited classics), the Ihara blueprint with the one missing part Z1 (Hilbert–Pólya in window coordinates), the measure and masses forced from pure lattice counting, and the corridor with the selection point at 0.53. Honest fence throughout: no claim of progress toward RH — roughly 700 machine-checked modules.",
  thumbnailUrl: [`${SITE_URL}${POSTER_SRC}`],
  uploadDate: "2026-08-03",
  duration: "PT2M39S",
  contentUrl: `${SITE_URL}${VIDEO_SRC}`,
  inLanguage: "de",
  isFamilyFriendly: true,
  transcript: TRANSCRIPT.map((s) => `${s.heading}\n${s.body}`).join("\n\n"),
  publisher: {
    "@type": "Organization",
    name: "TFPT Collaboration",
    url: SITE_URL,
  },
};

export function PrimeFrontVideo() {
  const ref = useRef<HTMLVideoElement>(null);
  const [started, setStarted] = useState(false);

  // Captions are burned into the video (open captions, German). Browsers
  // remember a user's prior "captions on" preference and would auto-show the
  // <track>, doubling the captions — so default the selectable track off;
  // users can still enable it from the native captions menu.
  const disableTracks = () => {
    const v = ref.current;
    if (!v) return;
    for (const t of Array.from(v.textTracks)) t.mode = "disabled";
  };
  useEffect(disableTracks, []);

  const play = () => {
    const v = ref.current;
    if (!v) return;
    setStarted(true);
    void v.play();
  };

  const seek = (t: number) => {
    const v = ref.current;
    if (!v) return;
    v.currentTime = t;
    setStarted(true);
    void v.play();
  };

  return (
    <section
      id="prime-front-video"
      className="relative scroll-mt-20 py-12 sm:py-16"
      aria-labelledby="prime-front-video-heading"
    >
      <script
        type="application/ld+json"
        dangerouslySetInnerHTML={{ __html: JSON.stringify(videoJsonLd) }}
      />
      <div className="mx-auto max-w-5xl px-4 sm:px-6 lg:px-8">
        <SectionHeader
          id="prime-front-video-heading"
          align="center"
          eyebrow="Watch first · auf Deutsch"
          title="Die Prime Front — in 2½ Minuten"
          description="The whole prime line as a short German-language film: RH in one image, the window matrix that turned out to be Suzuki's Weil operator, the calibrated detector, the two surface theorems, the Ihara blueprint with its one missing part — and the honest state: no proof, but the question has never been this small. German captions are burned in; the full transcript is below."
        />

        <figure className="mt-10">
          <div className="relative overflow-hidden border border-slate-700/50 bg-slate-950/60">
            <video
              ref={ref}
              className="aspect-video w-full"
              controls
              playsInline
              preload="none"
              poster={POSTER_SRC}
              aria-label="Die Prime Front — komplett erklärt (deutsch, mit Untertiteln)"
              onLoadedMetadata={disableTracks}
            >
              <source src={VIDEO_SRC} type="video/mp4" />
              {/* Captions are burned in (open captions); this selectable track
                  stays available for users who prefer browser-styled subtitles
                  and for indexing — not `default`, to avoid double captions. */}
              <track
                kind="subtitles"
                src={CAPTIONS_SRC}
                srcLang="de"
                label="Deutsch"
              />
              Your browser does not support the video tag. You can read the
              full transcript below.
            </video>

            {!started && (
              <button
                type="button"
                onClick={play}
                aria-label="Film abspielen — Die Prime Front, komplett erklärt"
                className="group absolute inset-0 flex items-center justify-center bg-slate-950/30 transition-colors hover:bg-slate-950/15 focus-visible:outline-none"
              >
                <span className="flex h-20 w-20 items-center justify-center rounded-full border border-sky-400/40 bg-sky-500/20 text-sky-100 ring-1 ring-white/10 transition-transform group-hover:scale-105 group-focus-visible:scale-105">
                  <Play size={34} className="ml-1" aria-hidden />
                </span>
              </button>
            )}
          </div>

          {/* chapter markers */}
          <nav
            aria-label="Video chapters"
            className="mt-4 flex flex-wrap justify-center gap-2"
          >
            {CHAPTERS.map((c) => (
              <button
                key={c.t}
                type="button"
                onClick={() => seek(c.t)}
                className={cn(
                  "inline-flex items-center gap-2 rounded-full border border-slate-700/50 bg-slate-900/60 px-3 py-1.5 text-xs font-medium text-slate-300 transition-colors hover:border-sky-400/40 hover:text-sky-200 focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-sky-400/40",
                )}
              >
                <span className="font-mono text-sky-300/80">{fmt(c.t)}</span>
                {c.label}
              </button>
            ))}
          </nav>

          <figcaption className="mt-3 text-center text-xs text-slate-500">
            German narration with German captions burned into the picture; a
            selectable subtitle track and the full transcript below are
            generated from the same source. Every number in the film traces to
            a named, machine-checked module — and the film claims no progress
            toward RH.
          </figcaption>
        </figure>

        {/* transcript — accessible + indexable */}
        <details className="group mt-8 rounded-2xl border border-slate-700/40 bg-slate-950/40 p-5 sm:p-6">
          <summary className="flex cursor-pointer list-none items-center justify-between text-sm font-semibold text-slate-200 transition-colors hover:text-white">
            <span>Transkript lesen (deutsch)</span>
            <span className="text-slate-500 transition-transform group-open:rotate-180">
              ⌄
            </span>
          </summary>
          <div className="prose-tfpt mt-5 space-y-5">
            {TRANSCRIPT.map((s) => (
              <div key={s.heading}>
                <h3 className="font-serif text-base font-semibold text-slate-100">
                  {s.heading}
                </h3>
                <p className="mt-1 text-sm leading-relaxed text-slate-300">
                  {s.body}
                </p>
              </div>
            ))}
          </div>
        </details>
      </div>
    </section>
  );
}
