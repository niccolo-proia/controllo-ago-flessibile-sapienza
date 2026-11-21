# Modellazione e controllo lineare della deflessione di un ago flessibile per procedure percutanee

Codice sorgente **MATLAB/Simulink** sviluppato nell'ambito della Tesi di Laurea Triennale in Ingegneria Informatica e Automatica (Sapienza Università di Roma, A.A. 2025/2026).

* **Autore:** Niccolò Proia
* **Relatore:** Prof.ssa Marilena Vendittelli

---

## Descrizione del Progetto

Il repository contiene i modelli matematici e gli schemi di simulazione per l'analisi e il controllo della deflessione planare di un ago flessibile in acciaio ($L = 0.15\text{ m}$) immerso in un tessuto molle modellato come fondazione elastica alla Winkler:

1. **Modellazione agli stati (FDM):** discretizzazione spaziale della trave di Eulero-Bernoulli in $N = 20$ segmenti (sistema LTI *stiff* di ordine 40).
2. **Analisi strutturale e sintesi continua:** analisi modale, mappe poli-zeri, confronto tra regolatore PID (puramente integrale) e LQR aumentato (41 stati), e confronto tra Osservatore di Luenberger e Filtro di Kalman.
3. **Simulatore Ibrido Sampled-Data (Simulink):** architettura continuo-discreta LQG con campionamento ecografico **Zero-Order Hold a 30 Hz** ($T_s = 0.033\text{ s}$), rumore bianco limitato in banda ($\text{Noise Power} = 10^{-9}$) e solutore implicito `ode15s`.
4. **Robustezza parametrica:** test di inseguimento su bersaglio a $5\text{ mm}$ con variazioni del $\pm 50\%$ sulla rigidezza tissutale $K_t$.

---

## Struttura del Repository

* `Modello_Ago_Flessibile.m` — Costruzione del modello FDM a 40 stati, analisi modale, mappe poli-zeri, sintesi PID/LQR e confronto Luenberger vs. Kalman nel tempo continuo.
* `tuning_LQR.m` — Ricalibrazione dei pesi LQG ($Q_{41,41} = 10$, $R = 50$) per l'architettura campionata a 30 Hz.
* `Modello_Tesi_Capitolo5_ZOH.slx` — Schema Simulink del sistema ibrido continuo-discreto.
* `test_rigido.m` — Aggiornamento della dinamica dell'impianto per tessuto rigido ($+50\% \, K_t$).
* `test_morbido.m` — Aggiornamento della dinamica dell'impianto per tessuto morbido ($-50\% \, K_t$).

---

## Requisiti e Ordine di Esecuzione

**Requisiti:** MATLAB e Simulink (testato su R2026b) con **Control System Toolbox**.

Per riprodurre correttamente le simulazioni senza errori di variabili mancanti nel Workspace, seguire questo ordine:

1. Eseguire `Modello_Ago_Flessibile.m` per generare le matrici del sistema (`A`, `Bu_base`, `C`), il Filtro di Kalman (`A_kalman`, `B_kalman`, `C_kalman`, `D_kalman`), la condizione iniziale `x0_stimato` e i grafici dei Capitoli 3 e 4.
2. Eseguire `tuning_LQR.m` per ricalibrare i guadagni `K_x` e `K_i` per il dominio campionato.
3. Aprire ed eseguire `Modello_Tesi_Capitolo5_ZOH.slx` (tempo di simulazione `10 s`, solver `ode15s`) per visualizzare la risposta nominale sullo Scope.
4. Eseguire `test_rigido.m` e rilanciare la simulazione Simulink per il test a $+50\% \, K_t$.
5. Eseguire `test_morbido.m` e rilanciare la simulazione Simulink per il test a $-50\% \, K_t$.

---

## Come Citare (BibTeX)

```bibtex
@misc{Proia2026NeedleSteering,
  author       = {Niccolò Proia},
  title        = {Modellazione e controllo lineare della deflessione di un ago flessibile per procedure percutanee},
  howpublished = {Tesi di Laurea Triennale, Sapienza Università di Roma},
  year         = {2026},
  url          = {https://github.com/niccolo-proia/controllo-ago-flessibile-sapienza}
}
