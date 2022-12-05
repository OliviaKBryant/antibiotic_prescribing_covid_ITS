# Antibiotic Prescribing in General Practice during COVID-19 and Beyond

Analysis code and data for an interrupted time-series analysis using the negative binomial model to assess the impact that the start of COVID-19 restrictions, March 2020, and the end of COVID-19 restrictions, July 2021, had on
antibiotic prescribing rates in general practice in England. Outcomes were measured as absolute number of items and the prescribing rate as the number of items prescribed per 100 appointments.

Analysis code is available in both R and SAS.

This analysis resulted in a letter published in The Lancet Infectious Diseases on [insert date] available at [inset DOI].
[Insert citation]


### Data Sources

Monthly prescription data was obtained from [OpenPrescribing](https://openprescribing.net/) from January 2018 to July 2022. 

Monthly number of GP appointments in England and appointment type was obtained from [NHS Digital](https://digital.nhs.uk/data-and-information/publications/statistical/appointments-in-general-practice) for the same months.

### Package Requirements

| Package     | Version    |
|------------|----------|
| tseries    | 0.10-52  |
| Epi        | 2.47     |
| MASS       | 7.3-58.1 |
| scales     | 1.1.1    |
| lubridate  | 1.9.0    |
| timechange | 0.1.1    |
| tsModel    | 0.6-1    |
| dplyr      | 1.0.1    |
| ggplot2    | 3.3.5    |
| patchwork  | 1.1.2    |
| ggtext     | 0.1.2    |
