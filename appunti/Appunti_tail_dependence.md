# Appunti di tesi - Tail Dependence

## Scopo del documento

Questo file raccoglie progressivamente:

- la mappa dei paper assegnati;
- l'ordine di lettura;
- definizioni e intuizioni;
- risultati importanti;
- collegamenti con Economia e Finanza;
- dubbi e domande da discutere con il relatore.

Stato attuale: **prima ricognizione comparativa**, basata su abstract, introduzioni, struttura e conclusioni/sezioni finali. Le formule e i risultati saranno verificati in dettaglio durante la lettura completa.

---

## 1. Idea generale dell'argomento

La tail dependence descrive quanto due o più variabili tendano a entrare **insieme** nelle rispettive code, cioè a registrare congiuntamente valori estremi.

In finanza la domanda tipica è: se un'attività registra una perdita eccezionalmente elevata, quanto è probabile che anche un'altra attività registri una perdita eccezionale? La correlazione ordinaria riassume la dipendenza complessiva, ma può non descrivere correttamente il comportamento congiunto durante crisi, default o forti movimenti di mercato.

Le copule sono centrali perché separano:

- le distribuzioni marginali delle singole variabili;
- la struttura di dipendenza che le unisce.

La tail dependence è quindi una proprietà della copula e non delle singole distribuzioni marginali.

### Prima distinzione essenziale

- **Coda inferiore**: eventi congiuntamente molto bassi. Per rendimenti finanziari, spesso corrisponde a forti perdite.
- **Coda superiore**: eventi congiuntamente molto alti. Il significato finanziario dipende dalla variabile usata: per rendimenti può indicare forti guadagni; per variabili di perdita, forti perdite.

Questa convenzione dovrà essere sempre dichiarata, perché il significato economico di upper e lower tail cambia se si modellano rendimenti oppure perdite.

---

## 2. Mappa comparativa dei tre testi

| Testo | Natura | Dimensione prevalente | Tema centrale | Ruolo probabile nel percorso |
|---|---|---:|---|---|
| Schmidt, *Tail Dependence* | Capitolo introduttivo e applicato | Bivariata | Definizione, calcolo e stima del coefficiente di tail dependence; applicazioni a rendimenti e VaR | Base concettuale e primo testo da leggere |
| Durante, Fernández-Sánchez, Pappadà, *Copulas, diagonals, and tail dependence* | Review teorico-applicata con nuovi strumenti | Soprattutto bivariata | Sezioni diagonali delle copule, tail concentration function, visualizzazione e clustering finanziario | Ponte tra nozioni di base, limiti del coefficiente e applicazioni |
| Hua e Joe, *Tail order and intermediate tail dependence of multivariate copulas* | Paper teorico avanzato | Multivariata | Tail order e tail order functions per distinguere dipendenza asintotica, intermedia e indipendenza | Approfondimento specialistico, da leggere per ultimo |

### Relazione logica tra i paper

I tre lavori non sembrano tre rami scollegati. Formano piuttosto una progressione:

1. **Schmidt** introduce il coefficiente classico di tail dependence, il suo calcolo e la sua stima.
2. **Durante et al.** mostra che un solo coefficiente asintotico può essere troppo povero e propone di osservare il comportamento lungo la diagonale anche a scala finita.
3. **Hua e Joe** affronta in modo più generale e tecnico un altro limite del coefficiente classico: copule con coefficiente nullo possono avere gradi molto diversi di associazione nelle code. Il tail order misura la velocità con cui decade la probabilità congiunta di coda.

Questa progressione suggerisce un possibile filo conduttore della tesi: **dal coefficiente classico a descrizioni più informative del comportamento congiunto nelle code**.

---

## 3. Ordine di lettura consigliato

### 1. Schmidt - *Tail Dependence*

È il testo più adatto per costruire il vocabolario di base. Parte dall'intuizione finanziaria, definisce i coefficienti di coda, passa alle copule e include calcolo, stima e applicazioni.

Obiettivi della lettura:

- capire intuitivamente e formalmente λL e λU;
- capire perché la correlazione non basta;
- vedere esempi di copule/distribuzioni tail-dependent e tail-independent;
- capire almeno a grandi linee come si stima la tail dependence;
- collegare la misura al VaR.

### 2. Durante et al. - *Copulas, diagonals, and tail dependence*

Va letto dopo avere acquisito le definizioni di base. Usa la diagonale della copula per descrivere le code e insiste sulla differenza tra comportamento asintotico e comportamento a soglie finite, importante nelle applicazioni reali.

Obiettivi della lettura:

- capire cos'è la sezione diagonale δC(t) = C(t,t);
- vedere come da essa derivano i coefficienti di coda;
- capire la tail concentration function;
- comprendere perché copule con lo stesso coefficiente limite possano comportarsi diversamente a livelli finiti;
- studiare le applicazioni grafiche e il clustering di serie finanziarie.

### 3. Hua e Joe - *Tail order and intermediate tail dependence of multivariate copulas*

È il paper più tecnico. Richiede familiarità con copule, limiti nelle code e, almeno intuitivamente, variazione regolare e funzioni lentamente variabili.

Obiettivi della lettura:

- capire il problema dell'indipendenza asintotica non uniforme;
- interpretare il tail order κ come velocità di decadimento della probabilità congiunta di coda;
- distinguere tail dependence, intermediate tail dependence e indipendenza;
- capire il ruolo delle tail order functions nelle diverse direzioni, non solo lungo la diagonale;
- seguire gli esempi sulle copule archimedee e sull'asimmetria tra coda inferiore e superiore.

---

## 4. Schede preliminari dei paper

### 4.1 Schmidt - *Tail Dependence*

#### Domanda centrale

Come si definisce, si calcola e si stima la dipendenza tra eventi estremi, e perché essa è rilevante per il rischio finanziario e il VaR?

#### Contenuto

Il capitolo:

1. introduce la tail dependence attraverso eventi estremi congiunti;
2. definisce i coefficienti superiore e inferiore;
3. li esprime tramite la copula;
4. li calcola per copule archimedee, distribuzioni ellittiche e altre famiglie;
5. presenta stimatori parametrici e non parametrici;
6. confronta empiricamente gli stimatori;
7. analizza rendimenti di attività e tassi di cambio;
8. svolge una simulazione sul Value at Risk.

#### Formule fondamentali iniziali

Per una coppia continua (X₁, X₂), il coefficiente di dipendenza nella coda superiore è:

```text
λU = lim [v → 1⁻] P(X₁ > F₁⁻¹(v) | X₂ > F₂⁻¹(v))
```

Si legge: quando la soglia percentile `v` si avvicina a 1, qual è la probabilità che X₁ superi il proprio quantile estremo, sapendo che X₂ ha già superato il suo?

Il coefficiente nella coda inferiore è

```text
λL = lim [v → 0⁺] P(X₁ ≤ F₁⁻¹(v) | X₂ ≤ F₂⁻¹(v))
```

In termini di copula C, le due formule diventano:

```text
                       1 - 2v + C(v,v)
λU = lim [v → 1⁻]     -----------------
                              1 - v

                              C(v,v)
λL = lim [v → 0⁺]             ------
                                 v
```

Queste formule anticipano il legame diretto con la diagonale studiata da Durante et al.

#### Esempio chiave

- La copula/distribuzione gaussiana con correlazione |ρ| < 1 è asintoticamente tail-independent.
- La distribuzione t bivariata può essere tail-dependent e il grado dipende da correlazione e gradi di libertà.

Quindi una correlazione elevata non implica necessariamente tail dependence positiva.

#### Livello di difficoltà

Medio e progressivo. È il testo più accessibile dei tre, anche se le parti sugli stimatori richiederanno attenzione.

---

### 4.2 Durante, Fernández-Sánchez e Pappadà - *Copulas, diagonals, and tail dependence*

#### Riferimento

Fuzzy Sets and Systems, 264 (2015), pp. 22-41.

#### Domanda centrale

Quanto si può conoscere o costruire della dipendenza di coda di una copula bivariata attraverso la sua sezione diagonale? Come descrivere le code in modo più ricco di un singolo coefficiente limite?

#### Struttura

1. Introduzione e motivazione nel risk management.
2. Proprietà e costruzioni di copule con diagonale assegnata.
3. Relazione tra diagonale e tail dependence.
4. Strumenti grafici:
   - visualizzazione dello spazio dei modelli di copula;
   - clustering di serie finanziarie in base al comportamento nelle code.
5. Conclusioni.

#### Contributo principale

Il paper combina una rassegna di risultati sulle copule con diagonale prescritta con la **tail concentration function**, pensata per descrivere la dipendenza a scala finita. Il messaggio centrale è che il limite asintotico λ da solo può nascondere differenze importanti tra modelli.

#### Collegamento diretto con Schmidt

Ponendo δC(t) = C(t,t):

- la coda inferiore dipende da δC(t) / t quando t → 0;
- la coda superiore dipende da [1 - 2t + δC(t)] / (1 - t) quando t → 1.

Sono precisamente le espressioni del coefficiente viste in Schmidt, ma Durante et al. osserva l'intera traiettoria prima del limite.

#### Messaggio applicativo

In campioni finanziari finiti non osserviamo mai veramente il limite. Occorre quindi considerare sia:

- la dipendenza asintotica;
- il comportamento pre-asintotico o a soglie finite.

La funzione proposta viene usata sia per confrontare famiglie di copule graficamente, sia per raggruppare serie finanziarie con comportamento di coda simile.

#### Livello di difficoltà

Medio-alto. La sezione 2 è più matematica; le sezioni 3 e 4 sono probabilmente le più centrali per una tesi in Economia e Finanza.

---

### 4.3 Hua e Joe - *Tail order and intermediate tail dependence of multivariate copulas*

#### Riferimento

Journal of Multivariate Analysis, 102 (2011), pp. 1454-1471.

#### Domanda centrale

Come descrivere la forza della dipendenza nelle code quando il coefficiente classico è zero, e come estendere l'analisi al caso multivariato e a direzioni diverse dalla diagonale?

#### Problema di partenza

Il coefficiente classico separa spesso i modelli in due classi:

- λ > 0: dipendenza asintotica;
- λ = 0: indipendenza asintotica.

Ma la seconda classe è molto ampia. Due modelli entrambi con λ = 0 possono produrre probabilità congiunte estreme che decadono a velocità molto diverse. Il tail order serve a misurare questa velocità.

#### Idea formale

Lungo la diagonale inferiore, il paper considera un comportamento del tipo

```text
C(u, ..., u) ∼ u^(κL) · ℓ(u),     per u → 0⁺
```

dove ℓ(u) è lentamente variabile. Un'espressione analoga vale per la coda superiore usando la survival copula.

Intuitivamente:

- κ = 1 corrisponde alla scala della tail dependence asintotica;
- 1 < κ < d descrive dipendenza intermedia;
- κ = d è la scala dell'indipendenza per una copula d-dimensionale;
- a parità di altre condizioni, un κ più piccolo indica un decadimento più lento e quindi eventi estremi congiunti relativamente più probabili.

Le tail order functions estendono l'analisi oltre la sola diagonale, studiando differenti direzioni verso l'angolo della copula.

#### Struttura

1. Motivazione tramite tail dependence, Ledford-Tawn e asimmetria di coda.
2. Definizioni e proprietà di tail orders e tail order functions.
3. Intermediate tail dependence per copule archimedee.
4. Copule costruite da mixture di copule max-infinitely-divisible.
5. Discussione e appendice con dimostrazioni.

#### Contributi specifici

- quadro integrato per tail dependence e intermediate tail dependence;
- proprietà multivariate del tail order;
- studio delle code di copule archimedee tramite trasformate di Laplace;
- nuova famiglia archimedea basata sull'inversa Gamma;
- analisi dell'asimmetria tra code inferiore e superiore;
- risultati per mixture di copule max-infinitely-divisible.

#### Livello di difficoltà

Alto. Richiede variazione regolare, funzioni lentamente variabili, omogeneità e trasformate di Laplace. Non è il testo con cui iniziare.

---

## 5. Concetti da padroneggiare prima dell'incontro

Priorità alta:

1. distribuzione congiunta e distribuzioni marginali;
2. quantili e funzioni inverse generalizzate;
3. probabilità condizionata;
4. copula e teorema di Sklar;
5. significato di upper e lower tail dependence;
6. coefficienti λU e λL;
7. differenza tra correlazione e dipendenza di coda;
8. dipendenza e indipendenza asintotica;
9. ruolo della diagonale C(u,u);
10. comportamento asintotico rispetto a comportamento a soglia finita.

Priorità successiva:

11. copule archimedee ed ellittiche;
12. survival copula;
13. tail concentration function;
14. intermediate tail dependence e tail order;
15. variazione regolare e funzioni lentamente variabili;
16. asimmetria tra coda inferiore e superiore.

---

## 6. Ipotesi iniziale sul possibile centro della tesi

Non è ancora possibile stabilire con certezza quale sarà il tema preciso della tesi. Tuttavia, la selezione dei tre testi suggerisce che il relatore possa essere interessato non soltanto al coefficiente classico, ma ai suoi limiti e a misure più ricche.

Possibili direzioni:

- confronto tra misure classiche e misure raffinate di dipendenza di coda;
- ruolo del comportamento a scala finita nella selezione di copule finanziarie;
- intermediate tail dependence in presenza di indipendenza asintotica;
- asimmetria delle code nei rendimenti finanziari;
- applicazioni a risk management, VaR, contagio, clustering o diversificazione di portafoglio;
- estensione dal caso bivariato al multivariato.

Queste sono ipotesi di lavoro, non conclusioni definitive.

---

## 7. Domande utili per il relatore

1. Il centro della tesi sarà teorico, empirico o una combinazione dei due?
2. Il paper di Hua e Joe rappresenta il nucleo teorico da sviluppare oppure un approfondimento dei limiti del coefficiente classico?
3. L'interesse principale è bivariato o multivariato?
4. È prevista un'applicazione a dati finanziari reali?
5. Quali famiglie di copule dovranno essere considerate?
6. Il relatore si aspetta stima non parametrica, stima parametrica o confronto tra entrambe?
7. Il comportamento a soglia finita della tail concentration function è una possibile direzione centrale?
8. Quale convenzione useremo: rendimenti oppure perdite? Di conseguenza, quale coda rappresenterà gli eventi avversi?

---

## 8. Piano di lavoro immediato

- [x] Identificare e classificare i tre testi.
- [x] Stabilire un primo ordine di lettura.
- [ ] Leggere Schmidt, sezioni 3.1-3.2: motivazione e definizioni.
- [ ] Costruire un glossario con esempi numerici/intuitivi.
- [ ] Leggere Schmidt, sezione 3.3: calcolo per famiglie di modelli.
- [ ] Leggere Schmidt, sezioni 3.4-3.7: stima e applicazioni.
- [ ] Passare a Durante et al., concentrandosi inizialmente sulle sezioni 1, 3, 4 e 5.
- [ ] Affrontare Hua e Joe dopo una breve preparazione sulla variazione regolare.

---

## 9. Nota metodologica sugli appunti

Durante le letture distingueremo sempre:

- **cosa afferma il paper**;
- **spiegazione intuitiva nostra**;
- **formula o risultato da ricordare**;
- **possibile rilevanza per la tesi**;
- **punto non ancora chiaro o da verificare**.

Questo eviterà di confondere il riassunto degli autori con le nostre interpretazioni preliminari.
