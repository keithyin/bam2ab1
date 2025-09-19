#[derive(Debug, Clone)]
pub struct PascalString {
    inner: Vec<u8>,
}

impl TryFrom<&str> for PascalString {
    type Error = String;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        if value.len() > 254 {
            return Err(format!("covert error. length>254"));
        }
        let mut inner = vec![0; value.len() + 1];
        inner[0] = value.len() as u8;
        inner[1..].copy_from_slice(value.as_bytes());
        Ok(Self { inner })
    }
}

impl PascalString {
    pub fn new(inner: Vec<u8>) -> Self {
        Self { inner }
    }

    pub fn bytes_len(&self) -> usize{
        self.inner.len()
    }

    pub fn bytes_raw(&self) -> Vec<u8> {
        self.inner.to_vec()
    }


    pub fn to_string(&self) -> anyhow::Result<String> {
        if self.inner.is_empty() {
            anyhow::bail!("not a valid pascal string");
        }
        let str_bytes = self.inner[0] as usize;
        if str_bytes > (self.inner.len() - 1) {
            anyhow::bail!("length exceed the actual value");
        }

        Ok(String::from_utf8(self.inner[1..].to_vec())?)
    }
}
