library(tidyverse)
library(nycflights13)

airlines
airports
planes
weather
flights %>% 
  count(year, month, day, flight) %>% 
  filter(n > 1)


#add surroagate key
flights %>% 
  arrange(year, month, day, sched_dep_time, carrier, flight) %>%
  mutate(flight_id = row_number()) %>%
  glimpse()


flights2 <- flights %>% 
  select(year:day, hour, origin, dest, tailnum, carrier)
flights2

flights2 %>%
  select(-origin, -dest) %>% 
  left_join(airlines, by = "carrier")
