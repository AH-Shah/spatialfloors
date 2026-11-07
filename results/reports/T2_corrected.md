| model | units | T mean | upper 95% | margin | locked = best | reading |
|---|---|---|---|---|---|---|
| F1|floor | 25 | +0.0043 | +0.0057 | 0.0195 | 0.19 | transfers |
| F1|scvi | 14 | +0.0021 | +0.0031 | 0.0195 | 0.20 | transfers |
| F1|nicheformer | 25 | +0.0001 | +0.0003 | 0.0195 | 0.92 | transfers |
| F2|floor | 10 | +0.0835 | +0.1356 | 0.0500 | 0.00 | inconclusive |
| F2|cellcharter | 5 | +0.0241 | +0.0766 | 0.0500 | 0.60 | inconclusive |
| F2|novae | 4 | -0.0000 | +0.0000 | 0.0500 | 1.00 | transfers |
| F3|floor | 22 | +0.0014 | +0.0020 | 0.0043 | 0.14 | transfers |
| F3|lightgbm | 22 | +0.0004 | +0.0008 | 0.0043 | 0.91 | transfers |
| F3|ncem | 11 | +0.0022 | +0.0049 | 0.0043 | 0.73 | inconclusive |
| F5|floor | 20 | +0.0025 | +0.0031 | 0.0255 | 0.05 | transfers |
| F5|scanvi | 11 | +0.0024 | +0.0043 | 0.0255 | 0.33 | transfers |
| F5|tangram | 11 | +0.0034 | +0.0070 | 0.0255 | 0.75 | transfers |

| pair | units | C mean | lower 95% | upper 95% | D same units | T learned | T floor | reading |
|---|---|---|---|---|---|---|---|---|
| F1|scvi | 14 | -0.0122 | -0.0162 | -0.0082 | -0.0093 | +0.0019 | +0.0047 | floor holds at best tuning |
| F1|nicheformer | 25 | -0.0421 | -0.0471 | -0.0371 | -0.0379 | +0.0001 | +0.0042 | floor holds at best tuning |
| F2|cellcharter | 5 | -0.2464 | -0.4262 | -0.0666 | -0.1242 | +0.0241 | +0.1463 | floor holds at best tuning |
| F2|novae | 4 | -0.1402 | -0.1724 | -0.1080 | -0.1147 | +0.0000 | +0.0255 | floor holds at best tuning |
| F3|lightgbm | 22 | -0.0166 | -0.0194 | -0.0137 | -0.0155 | +0.0004 | +0.0014 | floor holds at best tuning |
| F3|ncem | 11 | -0.0347 | -0.0429 | -0.0264 | -0.0352 | +0.0022 | +0.0017 | floor holds at best tuning |
| F5|scanvi | 11 | -0.0166 | -0.0297 | -0.0036 | -0.0170 | +0.0022 | +0.0021 | floor holds at best tuning |
| F5|tangram | 11 | -0.1183 | -0.1575 | -0.0791 | -0.1196 | +0.0045 | +0.0021 | floor holds at best tuning |

tier-1 report runs/p4/P5_corrected.json sha256 6d2b2b18339ae3a58e058d95a56a71ba2c28d5d19ec3e7ee933ab57f4abb36d3
