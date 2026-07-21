# Predicting the Next Mechanism State

FieldBridge separates two questions that are easy to confuse.

1. **Encoding:** Given an equation, which operational state does the frozen
   language assign to it?
2. **Prediction:** Given the current state and the first observed move, which
   move and operational state will the *next, unseen equation* occupy?

The word `destination` refers only to the second question. The destination is
the frozen-language assignment of a future equation that was withheld from the
model. It is not a second attempt to classify the equation already supplied.

The corpus-scale Hyperion evaluation reported the following complete-paper
held-out results:

| Withheld target | Top-1 accuracy |
|---|---:|
| Next move, `T2` | 57.4% |
| Future operator state, `Omega2` | 74.4% |
| Future substrate state, `Xi2` | 71.7% |
| Future completion state, `F2` | 73.4% |

The next-move model gained 0.210 bits per transition over
`P(T2 | T1)` (paper-cluster 95% interval 0.197 to 0.223 bits). The joint future
state gained 3.743 bits (95% interval 3.629 to 3.857 bits). The machine was
therefore predicting how a mechanism continued, rather than merely naming its
current state.

The frozen reference record is stored in
[`data/hyperion_future_state_reference.json`](../data/hyperion_future_state_reference.json).

## Public Reproduction Contract

`fieldbridge validate-continuation` applies the same distinction to a public
transition manifest. Each record describes two successive source-local moves:

```json
{
  "paper_id": "paper-001",
  "current_omega": "Omega05",
  "current_xi": "Xi03",
  "current_completion": "C+R",
  "first_move": "T05",
  "next_move": "T02",
  "destination_omega": "Omega05",
  "destination_xi": "Xi07",
  "destination_completion": "C+R+P"
}
```

Complete papers, not individual rows, are assigned to the evaluation fold.
The model receives `(current_omega, current_xi, current_completion, first_move)`
and predicts four withheld targets: `next_move`, `destination_omega`,
`destination_xi`, and `destination_completion`.

```bash
fieldbridge validate-continuation benchmark/transitions.jsonl \
  --min-evaluation-transitions 100 \
  --out-json build/future_state_validation.json \
  --out-md build/future_state_validation.md
```

The report gives top-1 and top-three accuracy and the code-length gain over a
first-move-only baseline for every target. This evaluates symbolic mechanism
continuation. It does not establish exact equation reconstruction or physical
validity.
