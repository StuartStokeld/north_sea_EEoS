# H1 results — plain-language guide to the statistics

**What this is for.** A non-technical walkthrough of the H1 headline numbers: what each statistic means, what the pattern shows, and how that should shape how the results are presented.

**Data in one line.** North Sea IBTS Q1 trawl hauls, 1985–2015, **12,069** hauls after joining survey length data to catch biomass. For each haul we ask: given species richness (S), abundance (N), and metabolic rate (E), does the Ecological Equation of State (EEoS) predict the haul’s biomass?

**Numbers below** match the current pipeline run (`outputs/h1_model_comparison.csv`).

---

## The question in plain English

EEoS is a theory that says community biomass should be predictable from three state variables — how many species, how many individuals, and how much metabolic “work” the community is doing.

**H1 asks:** on real North Sea survey hauls, does that prediction match what was actually caught?

Prior papers said yes, in a strong sense: EEoS biomass looked **about the same as a simple productivity map**, under a direct one-to-one comparison (no regression trained on the same dataset). H1 repeats that style of test for marine trawl data.

---

## How the test works (why “unfitted” matters)

Imagine plotting predicted biomass against observed biomass. A perfect theory would put every haul on the diagonal line “prediction = observation.”

Two rules make this a fair theory check:

1. **No free parameters tuned to this dataset.** We do not fit a slope or intercept to make EEoS look better on NS-IBTS. That would answer a different, easier question (“can we *calibrate* EEoS to these hauls?”).
2. **Same footing for the baseline.** The productivity comparison uses the same unfitted 1:1 line. Productivity’s scale constant (`m_min`) comes from the smallest fish in the haul — the same independent convention used to put EEoS predictions into grams — **not** from fitting to observed catch.

So the headline is not “which regression wins?” It is: **on a fair, parameter-free 1:1 map, does EEoS match catch biomass as well as a simple productivity proxy?**

---

## The three statistics you need

People often say “R²” as if there were one number. Here there are three, and they answer different questions.

### 1. Coefficient of determination — `log_r2` (the headline)

**Everyday meaning:** “How much better is this prediction than just guessing the average haul every time?”

- **1.0** = perfect  
- **0** = no better than predicting the mean  
- **Negative** = *worse* than predicting the mean (systematic mistakes beat a lazy constant guess)

It is computed on **log biomass**, so a prediction that is always ~4× too high is heavily punished — even if the ranking of hauls looks sensible.

| Model | log_r2 | Plain reading |
|-------|-------:|---------------|
| EEoS | **−0.22** | Worse than guessing the average. Absolute biomass prediction fails. |
| Productivity 1:1 (E × m_min) | **+0.74** | Strong unfitted map. Explains most log-scale variation. |

**Key takeaway:** under the same unfitted rules, a simple productivity proxy works; the full equation of state does not.

### 2. Correlation squared — `cor²` (diagnostic only)

**Everyday meaning:** “Do bigger predicted hauls tend to be bigger observed hauls?” Rank / pattern agreement, **ignoring** whether the scale is right.

| Model | cor² | Plain reading |
|-------|-----:|---------------|
| EEoS | **0.93** | Excellent ranking: EEoS knows which hauls are large vs small. |
| Productivity 1:1 | **0.90** | Also excellent ranking. |

**Why this matters.** EEoS can look “good” on cor² while failing on `log_r2`. That combination — tight cloud, wrong height — is the signature of a **systematic scale offset** (e.g. survey catchability), not random noise.

**Presentation rule:** never lead with cor² alone. It would hide a ~4× overprediction.

### 3. RMSE and median ratio (how wrong, and in which direction)

**log RMSE** — typical size of the log-scale error (lower is better).

| Model | log RMSE |
|-------|---------:|
| EEoS | **1.41** |
| Productivity 1:1 | **0.65** |

**Median B_pred / B_obs** — for a typical haul, how many times larger is the prediction than the catch?

| Model | Median pred/obs |
|-------|----------------:|
| EEoS | **4.1×** (overpredicts) |
| Productivity 1:1 | **0.67×** (slightly under) |

**Plain reading:** EEoS’s cloud sits systematically *above* the 1:1 line. The theory’s biomass is larger than what the survey caught, by about a factor of four at the median.

---

## The Harte criterion (one number for “does EEoS beat the baseline?”)

Prior work asked whether EEoS’s leftover error was less than half the baseline’s leftover error.

- **SS_res ratio (EEoS / productivity 1:1) = 4.6×**  
- Criterion (EEoS error &lt; 0.5 × productivity error): **not met**

**Plain reading:** EEoS leaves about **four and a half times** more unexplained log-scale variance than the simple unfitted productivity map. On this dataset, the equation of state does **not** earn the “equivalent to productivity” claim.

---

## Test 2 in one paragraph (relative pattern, not absolute biomass)

A second Harte figure asks whether the **ratio** E / B^(3/4) looks similar for predicted and observed biomass (relative structure, not absolute grams).

- Pearson cor² on raw ratios: **0.964** (all hauls) → **0.853** after dropping **2** extreme hauls  
- Harte’s published reference: **0.600**  
- Coefficient of determination on the same ratios is much lower (**~0.35**): pattern agreement without absolute-scale success

**Presentation rule:** report **both** all-hauls and trimmed; use trimmed as primary. Do not cite 0.964 alone. Relative-ratio agreement can look strong even when absolute biomass prediction fails — that is why Test 1 leads.

---

## What not to confuse with the headline

| Comparison | Role | Why it is secondary |
|------------|------|---------------------|
| Fitted ln(E) OLS (log_r2 ≈ **0.60**) | Extended correlative benchmark | Fits slope and intercept to *this* dataset — easier question; not like-for-like with unfitted theory checks |
| Uncalibrated E_raw 1:1 (log_r2 ≈ **0.62**) | Unit diagnostic only | Not the headline productivity baseline; superseded by E × m_min |

---

## Key takeaway for the analysis

1. **Absolute haul-level biomass prediction fails for EEoS** (`log_r2` negative; median ~4× overprediction).  
2. **A simple unfitted productivity map succeeds on the same footing** (`log_r2` ≈ +0.74).  
3. **High correlation does not rescue EEoS** — the theory ranks hauls well but misses the scale.  
4. **NS-IBTS does not behave like the terrestrial/arthropod systems** where prior work reported EEoS ≈ productivity under a direct 1:1 map.  
5. Residuals look **structured** (not random), which motivates H2/H3: do deviations track fishing pressure? — without claiming H1 already validated absolute prediction.

---

## Key takeaway for results presentation

**Lead with Test 1, unfitted, like-for-like:**

1. Show EEoS vs productivity 1:1 on the same 1:1 figure.  
2. Headline metric: **`log_r2`**, not cor².  
3. State the **4.1×** median overprediction and the **4.6×** SS_res ratio / Harte criterion not met.  
4. Only then: Test 2 ratio pattern (trimmed cor² primary), then fitted ln(E) as “different question / not comparable.”  
5. Spell out for any audience: *“Good ranking, wrong scale”* is the EEoS story; *“simple productivity proxy wins the unfitted comparison”* is the baseline story.

**One sentence you can say aloud:**  
On North Sea trawl hauls, EEoS does not predict catch biomass on an absolute scale, even though it ranks hauls well — and under the same fair 1:1 rules used in prior papers, a simple productivity map outperforms the full equation of state by a wide margin.
