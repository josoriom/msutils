use ionic::ion::{IonError, IonResult, Range};
use reqwest::{
    StatusCode,
    blocking::Client,
    header::{ACCEPT_RANGES, CONTENT_RANGE, RANGE},
};

pub struct UrlSource {
    client: Client,
    url: String,
}

impl UrlSource {
    pub fn new(url: String) -> IonResult<Self> {
        let client = Client::builder()
            .build()
            .map_err(|error| IonError::from(error.to_string()))?;
        Ok(Self { client, url })
    }

    pub fn read(&self, range: Range) -> IonResult<Vec<u8>> {
        let length = range.length;
        if length == 0 {
            return Ok(Vec::new());
        }

        let offset = range.offset;
        let last_byte = offset
            .checked_add(length)
            .and_then(|value| value.checked_sub(1))
            .ok_or_else(|| IonError::from("url source: range overflows"))?;
        let range = format!("bytes={offset}-{last_byte}");
        let response = self
            .client
            .get(&self.url)
            .header(RANGE, range)
            .send()
            .map_err(|error| IonError::from(error.to_string()))?;

        allow_range_response(response.status())?;
        allow_range_headers(&response, offset, last_byte)?;

        let bytes = response
            .bytes()
            .map_err(|error| IonError::from(error.to_string()))?
            .to_vec();
        allow_len(&bytes, length)?;
        Ok(bytes)
    }
}

fn allow_range_response(status: StatusCode) -> IonResult<()> {
    if status == StatusCode::PARTIAL_CONTENT {
        return Ok(());
    }
    Err(IonError::from(format!(
        "url source: expected 206 Partial Content, got {status}"
    )))
}

fn allow_range_headers(
    response: &reqwest::blocking::Response,
    offset: u64,
    last_byte: u64,
) -> IonResult<()> {
    if let Some(value) = response.headers().get(ACCEPT_RANGES)
        && value
            .to_str()
            .unwrap_or_default()
            .eq_ignore_ascii_case("none")
    {
        return Err(IonError::from("url source: server rejects range reads"));
    }

    let expected = format!("bytes {offset}-{last_byte}/");
    if response
        .headers()
        .get(CONTENT_RANGE)
        .and_then(|value| value.to_str().ok())
        .unwrap_or_default()
        .starts_with(&expected)
    {
        return Ok(());
    }
    Err(IonError::from("url source: bad content-range header"))
}

fn allow_len(bytes: &[u8], length: u64) -> IonResult<()> {
    let expected_len =
        usize::try_from(length).map_err(|_| IonError::from("url source: length out of range"))?;
    if bytes.len() == expected_len {
        return Ok(());
    }
    Err(IonError::from(format!(
        "url source: expected {expected_len} bytes, got {}",
        bytes.len()
    )))
}
